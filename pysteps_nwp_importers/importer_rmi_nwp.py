"""
pysteps_nwp_importers.importer_rmi_nwp
====================

Module to import the RMI NWP forecasts. The output of this method is a numpy
array containing forecast rainfall fields in mm/timestep and the metadata as
dictionary.

The metadata contains the following key-value pairs:

.. tabularcolumns:: |p{2cm}|L|

+------------------+----------------------------------------------------------+
|       Key        |                Value                                     |
+==================+==========================================================+
|   projection     | PROJ.4-compatible projection definition                  |
+------------------+----------------------------------------------------------+
|   x1             | x-coordinate of the lower-left corner of the data raster |
+------------------+----------------------------------------------------------+
|   y1             | y-coordinate of the lower-left corner of the data raster |
+------------------+----------------------------------------------------------+
|   x2             | x-coordinate of the upper-right corner of the data raster|
+------------------+----------------------------------------------------------+
|   y2             | y-coordinate of the upper-right corner of the data raster|
+------------------+----------------------------------------------------------+
|   xpixelsize     | grid resolution in x-direction                           |
+------------------+----------------------------------------------------------+
|   ypixelsize     | grid resolution in y-direction                           |
+------------------+----------------------------------------------------------+
|   cartesian_unit | the physical unit of the cartesian x- and y-coordinates: |
|                  | e.g. 'm' or 'km'                                         |
+------------------+----------------------------------------------------------+
|   yorigin        | a string specifying the location of the first element in |
|                  | the data raster w.r.t. y-axis:                           |
|                  | 'upper' = upper border                                   |
|                  | 'lower' = lower border                                   |
+------------------+----------------------------------------------------------+
|   institution    | name of the institution who provides the data            |
+------------------+----------------------------------------------------------+
|   unit           | the physical unit of the data: 'mm/h', 'mm' or 'dBZ'     |
+------------------+----------------------------------------------------------+
|   transform      | the transformation of the data: None, 'dB', 'Box-Cox' or |
|                  | others                                                   |
+------------------+----------------------------------------------------------+
|   accutime       | the accumulation time in minutes of the data, float      |
+------------------+----------------------------------------------------------+
|   threshold      | the rain/no rain threshold with the same unit,           |
|                  | transformation and accutime of the data.                 |
+------------------+----------------------------------------------------------+
|   zerovalue      | the value assigned to the no rain pixels with the same   |
|                  | unit, transformation and accutime of the data.           |
+------------------+----------------------------------------------------------+

"""

import numpy as np

import xarray as xr

from pysteps_nwp_importers.exceptions import MissingOptionalDependency

try:
    import netCDF4

    NETCDF4_IMPORTED = True
except ImportError:
    NETCDF4_IMPORTED = False

try:
    import dask

    DASK_IMPORTED = True
except ImportError:
    DASK_IMPORTED = False


def import_rmi_nwp(filename, **kwargs):
    """Import a NetCDF with NWP rainfall forecasts from RMI using xarray.

    Parameters
    ----------
    filename: str
        Name of the file to import.

    {extra_kwargs_doc}

    Returns
    -------
    precipitation : array-like, float32
        Precipitation field in mm/h. The dimensions are [time, rows, cols].
    quality : 2D array or None
        If no quality information is available, set to None.
    metadata : dict
        Associated metadata (pixel sizes, map projections, etc.).
    """

    if not NETCDF4_IMPORTED:
        raise MissingOptionalDependency(
            "netCDF4 package is required to import RMI NWP rainfall "
            "products but it is not installed"
        )

    if not DASK_IMPORTED:
        raise MissingOptionalDependency(
            "dask package is required to import RMI NWP regridded rainfall "
            "products but it is not installed"
        )

    ds = _import_rmi_nwp_data_xr(filename, **kwargs)
    metadata = _import_rmi_nwp_geodata_xr(ds, **kwargs)

    # rename varname_time (def: time) to t
    varname_time = kwargs.get("varname_time", "time")
    ds = ds.rename({varname_time: "t"})
    varname_time = "t"
    times_nwp = ds["t"].values
    metadata["time_stamps"] = times_nwp

    # if data variable is named accum_prcp
    # it is assumed that NWP rainfall data is accumulated
    # so it needs to be disagregated by time step
    varname = kwargs.get("varname", "precipitation")
    if varname == "accum_prcp":
        print("Rainfall values are accumulated. Disaggregating by time step")
        accum_prcp = ds[varname]
        precipitation = accum_prcp - accum_prcp.shift({varname_time: 1})
        precipitation = precipitation.dropna(varname_time, how="all")
        # update/copy attributes
        precipitation.name = "precipitation"
        # copy attributes
        precipitation.attrs.update({**accum_prcp.attrs})
        # update attributes
        precipitation.attrs.update({"standard_name": "precipitation_amount"})
    else:
        precipitation = ds[varname]

    quality = None

    return precipitation.values, quality, metadata


def _import_rmi_nwp_data_xr(filename, **kwargs):
    return xr.open_dataset(filename)


def _import_rmi_nwp_geodata_xr(
    ds_in,
    **kwargs,
):
    varname = kwargs.get("varname", "precipitation")
    varname_time = kwargs.get("varname_time", "time")
    projdef = None
    if "proj4string" in ds_in.attrs:
        projdef = ds_in.proj4string

    # get the accumulation period
    time = ds_in[varname_time]
    # shift the array to calculate the time step
    delta_time = time - time.shift({varname_time: 1})
    # assuming first valid delta_time is representative of all time steps
    time_step = delta_time[1]
    time_step = time_step.values.astype("timedelta64[m]").astype(int)

    # get the units of precipitation
    units = None
    if "units" in ds_in[varname].attrs:
        units = ds_in[varname].units
        if units in ("kg m-2", "mm"):
            units = "mm"

    # get spatial boundaries and pixelsize
    # move to meters if coordinates in kilometers
    if "units" in ds_in.x.attrs:
        if ds_in.x.units == "km":
            ds_in["x"] = ds_in.x * 1000.0
            ds_in.x.attrs.update({"units": "m"})
            ds_in["y"] = ds_in.y * 1000.0
            ds_in.y.attrs.update({"units": "m"})

    xmin = ds_in.x.min().values
    xmax = ds_in.x.max().values
    ymin = ds_in.y.min().values
    ymax = ds_in.y.max().values
    xpixelsize = abs(ds_in.x[1] - ds_in.x[0])
    ypixelsize = abs(ds_in.y[1] - ds_in.y[0])

    # Add metadata needed by pySTEPS as attrs in rainfall variable
    da_rainfall = ds_in[varname].isel({varname_time: 0})

    # Fill the metadata dictionary
    metadata = dict(
        xpixelsize=xpixelsize.values,
        ypixelsize=ypixelsize.values,
        cartesian_unit=ds_in.x.units,
        unit=units,
        transform=None,
        zerovalue=np.nanmin(da_rainfall),
        institution="Royal Meteorological Institute of Belgium",
        projection=projdef,
        yorigin="upper",
        threshold=_get_threshold_value(da_rainfall.values),
        x1=xmin,
        x2=xmax,
        y1=ymin,
        y2=ymax,
        accutime=time_step,
    )

    return metadata


def _get_threshold_value(precip):
    """
    Get the rain/no rain threshold with the same unit, transformation and
    accutime of the data.
    If all the values are NaNs, the returned value is `np.nan`.
    Otherwise, np.min(precip[precip > precip.min()]) is returned.

    Returns
    -------
    threshold: float
    """
    valid_mask = np.isfinite(precip)
    if valid_mask.any():
        _precip = precip[valid_mask]
        min_precip = _precip.min()
        above_min_mask = _precip > min_precip
        if above_min_mask.any():
            return np.min(_precip[above_min_mask])
        else:
            return min_precip
    else:
        return np.nan
