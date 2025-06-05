"""
pysteps_nwp_importers.importer_dwd_nwp
====================

Module to import the DWD ICON-RUC NWP forecasts. The output of this method
is a xarray containing the desired precipitation related variable per
timestep and the metadata as dictionary.

In case of rain rate, acc. prec., etc., the data is on an unstructured
triangular grid and the metadata contain the following key-value pairs:

.. tabularcolumns:: |p{2cm}|L|

+------------------+----------------------------------------------------------+
|       Key        |                Value                                     |
+==================+==========================================================+
|   projection     | PROJ.4-compatible projection definition                  |
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

In case of EMVORADO synthetic reflectivities, the data is on a regular grid
and the metadata contain additionally the following key-value pairs:

.. tabularcolumns:: |p{2cm}|L|

+------------------+----------------------------------------------------------+
|       Key        |                Value                                     |
+==================+==========================================================+
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

"""

import datetime
import pyproj
import numpy as np
import xarray as xr
import netCDF4 as nc

try:
    import pygrib

    PYGRIB_IMPORTED = True
except ImportError:
    PYGRIB_IMPORTED = False

from pysteps_nwp_importers.exceptions import MissingOptionalDependency


def import_dwd_nwp(filename, **kwargs):
    """Import a GRIB with ICON-RUC NWP precipitation forecasts from DWD
    using pygrib and xarray.

    Parameters
    ----------
    filename: str
        Name of the file to import.

    {extra_kwargs_doc}

    Returns
    -------
    precipitation : array-like, float32
        Either syn. reflectivity in dBZ with dimensions [time, rows, cols]
        or rain rate or acc. prec in mm/h and mm, resptively. The
        dimensions are [time, cols].
    quality : 2D array or None
        If no quality information is available, set to None.
    metadata : dict
        Associated metadata (pixel sizes, map projections, etc.).
    """

    if not PYGRIB_IMPORTED:
        raise MissingOptionalDependency(
            "pygrib package is required to import DWD NWP precipitation "
            "forecasts but it is not installed"
        )

    # try to open file with pygrib
    try:
        grib_file = pygrib.open(filename)
    except Exception as e:
        raise IOError("File could not be opened, because of: " f"{e}")

    # read grib file
    grib_msgs = grib_file.read()
    grib_file.close()

    # check whether the file contains any processible precipitation variable
    varname = kwargs.get("varname", "DBZCMP_SIM")
    data_varnames = np.unique([grib_msg["shortName"] for grib_msg in grib_msgs])

    if not varname in data_varnames:
        raise IOError("File does not contain the desired precipitation variable")

    assert (
        len(data_varnames) == 1
    ), "File should contain only one precipitation variable."

    reference_time = np.unique([grib_msg.analDate for grib_msg in grib_msgs])
    assert len(reference_time) == 1, "File should contain only one forecast init."

    if (
        grib_msgs[0].has_key("typeOfStatisticalProcessing")
        and grib_msgs[0]["typeOfStatisticalProcessing"] == 1
    ):
        laccum = True
        valid_times = [_valid_time_helper(grib_msg) for grib_msg in grib_msgs]
    else:
        laccum = False
        valid_times = [grib_msg.validDate for grib_msg in grib_msgs]

    ens_no = np.unique(
        [
            (
                str(grib_msg["perturbationNumber"])
                if grib_msg.has_key("perturbationNumber")
                else "0"
            )
            for grib_msg in grib_msgs
        ]
    )

    assert len(valid_times) / len(ens_no) == len(
        np.unique(valid_times)
    ), "No. of time stamps not equal between ensemble members."
    valid_times = np.unique(valid_times)

    time_steps = [
        int((valid_time - reference_time[0]).total_seconds() / 60)
        for valid_time in valid_times
    ]
    temp_res = np.diff(time_steps)[0]

    # get metadata and create DataArray
    da_prec, metadata = _import_dwd_nwp_geodata(
        grib_msgs[0], valid_times, ens_no, **kwargs
    )

    # fill DataArray
    for grib_msg in grib_msgs:
        valid_time = _valid_time_helper(grib_msg) if laccum else grib_msg.validDate
        ens_no = (
            str(grib_msg["perturbationNumber"])
            if grib_msg.has_key("perturbationNumber")
            else "0"
        )
        da_prec.loc[dict(time=valid_time, ens_no=ens_no)] = grib_msg["values"]

    metadata["time_stamps"] = valid_times
    # unfortunately, threshold values for synthetic reflectivities are hard coded here
    metadata["zerovalue"] = (
        np.nanmin(da_prec) if varname != "DBZCMP_SIM" else -2.0
    )  # grib_msgs[0]["referenceValue"]
    metadata["threshold"] = (
        _get_threshold_value(da_prec.to_numpy()) if varname != "DBZCMP_SIM" else 2.0
    )
    metadata["accutime"] = temp_res if laccum else None

    quality = None

    return da_prec, quality, metadata


def _import_dwd_nwp_geodata(grib_msg, valid_times, ens_no, **kwargs):

    units = None
    if "units" in grib_msg.keys():
        units = grib_msg["units"]
        if units in ("kg m-2", "mm"):
            units = "mm"

    # unfortunately projection is hard coded here since projparams does not support unstructured grids
    proj_params = grib_msg.projparams
    if proj_params is None:
        proj_def = (
            "+a=6370040.0 +b=6370040.0 +proj=stere +lat_ts=60.0 +lat_0=90.0 +lon_0=10.0"
        )
    else:
        proj_def = " ".join([f"+{key}={value} " for key, value in proj_params.items()])

    if grib_msg.has_key("Nj"):
        dims = ["time", "ens_no", "south_north", "west_east"]
        _, lat, lon = grib_msg.data()
        proj = pyproj.Proj(proj_params)
        x, y = proj(lon, lat)
        x = x[0]
        y = y[:, 0]
        coords = {
            "time": valid_times,
            "ens_no": ens_no,
            "west_east": np.arange(0, len(x), 1),
            "south_north": np.arange(0, len(y), 1),
            "x": ("west_east", np.arange(0, len(x), 1), {"units": "1"}),
            "y": ("south_north", np.arange(0, len(y), 1), {"units": "1"}),
            "projection_x_coordinate": ("west_east", x, {"units": "m"}),
            "projection_y_coordinate": ("south_north", y, {"units": "m"}),
            "longitude": (
                ["south_north", "west_east"],
                lon,
                {"units": "degrees_north"},
            ),
            "latitude": (["south_north", "west_east"], lat, {"units": "degrees_east"}),
        }
        da_prec = xr.DataArray(
            data=np.full(
                (len(valid_times), len(ens_no), grib_msg["Nj"], grib_msg["Ni"]), np.nan
            ),
            dims=dims,
            coords=coords,
        )

        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        xpixelsize = abs(x[1] - x[0])
        ypixelsize = abs(y[1] - y[0])

        metadata = dict(
            xpixelsize=xpixelsize,
            ypixelsize=ypixelsize,
            cartesian_unit="m",
            yorigin="lower",  # thats defined in section 3
            x1=xmin,
            x2=xmax,
            y1=ymin,
            y2=ymax,
        )

    else:
        grid_file_path = kwargs.get("grid_file_path", "")
        clon, clat = _read_grid_file(grid_file_path)
        dims = ["time", "ens_no", "cell_no"]
        coords = {
            "time": valid_times,
            "ens_no": ens_no,
            "cell_no": np.arange(0, len(clat), 1),
            "longitude": (["cell_no"], clon, {"units": "degrees_north"}),
            "latitude": (["cell_no"], clat, {"units": "degrees_east"}),
        }
        da_prec = xr.DataArray(
            data=np.full(
                (len(valid_times), len(ens_no), grib_msg["numberOfDataPoints"]), np.nan
            ),
            dims=dims,
            coords=coords,
        )
        metadata = dict()

    metadata["projection"] = proj_def
    metadata["unit"] = units
    metadata["transform"] = None
    metadata["institution"] = grib_msg["centre"]

    return da_prec, metadata


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


def _valid_time_helper(grib_msg):
    """
    Pygrib does not set the valid time correctly in case of total
    precipitation. Therefore, this helper function is necessary

    Returns
    -------
    validTime: datetime.datetime
    """
    validTime = datetime.datetime(
        grib_msg["yearOfEndOfOverallTimeInterval"],
        grib_msg["monthOfEndOfOverallTimeInterval"],
        grib_msg["dayOfEndOfOverallTimeInterval"],
        grib_msg["hourOfEndOfOverallTimeInterval"],
        grib_msg["minuteOfEndOfOverallTimeInterval"],
        grib_msg["secondOfEndOfOverallTimeInterval"],
    )
    return validTime


def _read_grid_file(grid_file_path):
    """
    Reads the coordinates of the unstructured triangular grid
    from the grid file

    Returns
    -------
    lat,lon: float
    """
    ds = nc.Dataset(grid_file_path)

    # center lon and lat of the triangles
    # -1 since these are fortran-like indices
    clon = ds.variables["clon"][:] - 1
    clat = ds.variables["clat"][:] - 1

    return clon, clat
