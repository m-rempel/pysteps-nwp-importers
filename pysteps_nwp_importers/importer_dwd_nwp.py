"""
pysteps_nwp_importers.importer_dwd_nwp
====================

Module to import the DWD ICON-RUC NWP forecasts. The output of this method
is an xarray DataArray containing the desired precipitation related variable per
timestep as well as the metadata as a dictionary.

The description of the metadata has to divide in two categories since
rainfall related variables are available only on an unstructured triangular
grid whereas EMVORADO-simulated reflectivities are available on a rotated
lat-lon grid.

Thus, for rain rate, acc. precipitation, etc., the metadata contain the
following key-value pairs:

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

In case of EMVORADO synthetic reflectivities, the metadata contain additionally
following entries:

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

A comprehensive documentation of the data and the NWP model itself can be test
found at: https://www.dwd.de/SharedDocs/downloads/DE/modelldokumentationen/nwv/
    icon_d2/icon_d2_dbbeschr_aktuell.pdf?view=nasPublication&nn=346850

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
    """
    Import a GRIB file containing DWD ICON-D2-RUC NWP precipitation forecasts using
    pygrib and xarray.

    Parameters
    ----------
    filename: str
        Name of the file to import.

    Other Parameters
    ----------------
    varname: str
        GRIB short name of the desired variable (depends on the installed GRIB tables):
        - DBZCMP_SIM = EMVORADO synthetic reflectivities
        - PR_GSP = grid scale precipitation rate
        - TOT_PREC = accumulated precipitation since forecast initialization
    grid_file_path: str
        Path to the forecast data associated grid file

    Returns
    -------
    tuple
        A tuple containing the following:
        - precipitation : xarray.DataArray, float32
            Either syn. reflectivity in dBZ with dimensions [time, rows, cols]
            or rain rate or acc. prec in mm/h and mm, resptively. Dimensions of the
            latter are [time, cols].
        - quality : np.ndarray or None
            If there's no quality information, it is set to None.
        - metadata : dict
            Associated metadata (pixel sizes, map projections, etc.; see doc of
            pysteps_nwp_importers.importer_dwd_nwp).
    """

    if not PYGRIB_IMPORTED:
        raise MissingOptionalDependency(
            "pygrib package is required to import DWD NWP precipitation "
            "forecasts but it is not installed"
        )

    # Try to open file with pygrib
    try:
        grib_file = pygrib.open(filename)
    except Exception as e:
        raise IOError("File could not be opened, because of: " f"{e}")

    # Read complete grib file
    grib_msgs = grib_file.read()
    grib_file.close()

    # Check whether the file contains any processible precipitation variable
    varname = kwargs.get("varname", "DBZCMP_SIM")
    data_varnames = np.unique([grib_msg["shortName"] for grib_msg in grib_msgs])

    # Raise error if file doesn't contain the desired variable
    if not varname in data_varnames:
        raise IOError("File does not contain the desired precipitation variable")

    # Raise error if file contains more than one variable
    assert (
        len(data_varnames) == 1
    ), "File should contain only one precipitation variable."

    # Check the initialization time
    reference_time = np.unique([grib_msg.analDate for grib_msg in grib_msgs])
    # Raise error if file contains more than one forecast init
    assert len(reference_time) == 1, "File should contain only one forecast init."

    # Check whether the file contains an accumulated variable
    if (
        grib_msgs[0].has_key("typeOfStatisticalProcessing")
        and grib_msgs[0]["typeOfStatisticalProcessing"] == 1
    ):
        laccum = True
        valid_times = [_valid_time_helper(grib_msg) for grib_msg in grib_msgs]
    else:
        laccum = False
        valid_times = [grib_msg.validDate for grib_msg in grib_msgs]

    # Get all ensemble members that the dataset contains
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

    # Check whether each ensemble member has the same number of forecast times
    # Currently, only the number is compared and not the individual forecast times
    # themselves.
    assert len(valid_times) / len(ens_no) == len(
        np.unique(valid_times)
    ), "No. of time stamps not equal between ensemble members."
    valid_times = np.unique(valid_times)

    # Get time steps and temporal resolution
    time_steps = [
        int((valid_time - reference_time[0]).total_seconds() / 60)
        for valid_time in valid_times
    ]
    temp_res = np.diff(time_steps)[0]

    # Get metadata and initialize DataArray
    da_prec, metadata = _import_dwd_nwp_geodata(
        grib_msgs[0], valid_times, ens_no, **kwargs
    )

    # Fill DataArray with values of the grib messages by ensemble member and forecast
    # time
    for grib_msg in grib_msgs:
        valid_time = _valid_time_helper(grib_msg) if laccum else grib_msg.validDate
        ens_no = (
            str(grib_msg["perturbationNumber"])
            if grib_msg.has_key("perturbationNumber")
            else "0"
        )
        da_prec.loc[dict(time=valid_time, ens_no=ens_no)] = grib_msg["values"]

    # Set forecast times additionally in the metadate
    metadata["time_stamps"] = da_prec["time"].values
    # Set threshold and zerovalue of the dataset
    # For synthetic reflectivities no data and noch echo is set in the GRIB to -999 and
    # -99, respectively. To handle these in a more comfortable manner, it is set here
    # to -2 and 2. However, currently hard coded...
    metadata["zerovalue"] = np.nanmin(da_prec) if varname != "DBZCMP_SIM" else -2.0
    metadata["threshold"] = (
        _get_threshold_value(da_prec.to_numpy()) if varname != "DBZCMP_SIM" else 2.0
    )
    # Set the accumulation time equal to time steps if it is a variable that is
    # accumulated since forecast initialization. Set it to the temporal resolution if it is a variable accumulated within two time steps. Otherwise set it to None.
    if laccum and "TOT_" in varname:
        metadata["accutime"] = time_steps
    elif laccum:
        metadata["accutime"] = temp_res
    else:
        metadata["accutime"] = None

    # There is no information about quality. Therefore, for consistence reasons,
    # quality is set to None.
    quality = None

    return da_prec, quality, metadata


def _import_dwd_nwp_geodata(grib_msg, valid_times, ens_no, **kwargs):
    """
    Get all necessary metadata for further processing from data the GRIB file contains.

    Parameters
    ----------
    grib_msg: Grib message object
        Single grib message from which metadata is to be extracted
    valid_times: np.ndarray
        1D array containing the list of included forecast times
    ens_no: np.ndarray
        1D array containing the list of included ensemble members

    Other parameters
    ----------------
    grid_file_path: str
        Path to the forecast data associated grid file

    Returns
    -------
    tuple
        A tuple containing the following:
        - da_prec : xarray.DataArray, float32
            With geodata initialized DataArray
        - metadata : dict
            Associated metadata (pixel sizes, map projections, etc.; see doc of
            pysteps_nwp_importers.importer_dwd_nwp).
    """

    # Set the unit of the desired variable
    units = None
    if "units" in grib_msg.keys():
        units = grib_msg["units"]
        if units in ("kg m-2", "mm"):
            units = "mm"

    # For the rotated lat/lon grid the projection definition is extracted from
    # pygrib.proj_params. For the unstructured grid, it is unfortunately hard coded,
    # since proj_params does not support this kind of grids.
    proj_params = grib_msg.projparams
    if proj_params is None:
        proj_def = (
            "+a=6378137.0 +b=6356752.0 +proj=stere +lat_ts=60.0 +lat_0=90.0 +lon_0=10.0"
        )
    else:
        proj_def = " ".join([f"+{key}={value} " for key, value in proj_params.items()])

    # Initialize DataArray and get metadata for a rotated lat/lon grid
    if grib_msg.has_key("Nj"):
        # Set dimensions for DataArray
        dims = ["time", "ens_no", "south_north", "west_east"]
        # Get geographical coordinates from grib message
        _, lat, lon = grib_msg.data()
        # Get projection and corresponding carthesian coordinates
        proj = pyproj.Proj(proj_params)
        x, y = proj(lon, lat)
        x = x[0]
        y = y[:, 0]
        # Fill coordinate dictionary for DataArray
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
        # Initiliaze DataArray
        da_prec = xr.DataArray(
            data=np.full(
                (len(valid_times), len(ens_no), grib_msg["Nj"], grib_msg["Ni"]), np.nan
            ),
            dims=dims,
            coords=coords,
        )

        # Get corners of the carthesian grid and set orientation of y-axis
        xmin = x.min()
        xmax = x.max()
        ymin = y.min()
        ymax = y.max()
        xpixelsize = abs(x[1] - x[0])
        ypixelsize = abs(y[1] - y[0])
        yorigin = "upper" if grib_msg["orientationOfTheGrid"] < 0 else "lower"

        # Fill the metadata dictionary with horizontal resolution, length unit,
        # orientation of y-axis as well as corner coordinates
        metadata = dict(
            xpixelsize=xpixelsize,
            ypixelsize=ypixelsize,
            cartesian_unit="m",
            yorigin=yorigin,
            x1=xmin,
            x2=xmax,
            y1=ymin,
            y2=ymax,
        )

    # Initialize DataArray and get metadata for an unstructured triangular grid
    else:
        # Set dimensions for DataArray
        dims = ["time", "ens_no", "cell_no"]
        # Get geographical coordinates of triangle centers from grid file
        grid_file_path = kwargs.get("grid_file_path", "")
        clon, clat = _read_grid_file(grid_file_path)
        # Fill coordinate dictionary for DataArray
        coords = {
            "time": valid_times,
            "ens_no": ens_no,
            "cell_no": np.arange(0, len(clat), 1),
            "longitude": (["cell_no"], clon, {"units": "degrees_north"}),
            "latitude": (["cell_no"], clat, {"units": "degrees_east"}),
        }
        # Initiliaze DataArray
        da_prec = xr.DataArray(
            data=np.full(
                (len(valid_times), len(ens_no), grib_msg["numberOfDataPoints"]), np.nan
            ),
            dims=dims,
            coords=coords,
        )
        metadata = dict()

    # Fill metadata dictionary with projection definition, unit of the precipitation
    # variable, transformation as well as the originating institution
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

    Parameters
    ----------
    precip: np.ndarray
        2D array containing the precipitation field from which the threshold value
        is to be extracted

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

    Parameters
    ----------
    grib_msg: Grib message object
        GRIB message of the file to read

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

    Parameters
    ----------
    grid_file_path: str
        Path to the forecast data associated grid file

    Returns
    -------
    lat,lon: float
    """
    # Open grid file
    ds = nc.Dataset(grid_file_path)

    # Get lon/lat coordinates of triagle centers of the grid
    clon = np.rad2deg(ds.variables["clon"][:])
    clat = np.rad2deg(ds.variables["clat"][:])

    return clon, clat
