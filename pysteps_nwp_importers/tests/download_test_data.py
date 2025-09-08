import os
from pathlib import Path
from urllib import request
from urllib.parse import urljoin

DATA_DIR = Path(__file__).parent / "../tests/data"


def download_test_data():
    root_url = "https://github.com/pySTEPS/pysteps-data/raw/master/"

    files_to_download = (
        "nwp/bom/2020/10/31/20201031_0000_regrid_short.nc",
        "nwp/knmi/2018/09/05/20180905_0600_Pforecast_Harmonie.nc",
        "nwp/rmi/2021/07/04/ao13_2021070412_native_5min.nc",
        "nwp/dwd/2025/06/04/20250604_1600_PR_GSP_060_120.grib2",
        "aux/grid_files/dwd/icon/R19B07/icon_grid_0047_R19B07_L.nc",
    )

    for _file in files_to_download:
        file_url = urljoin(root_url, _file)
        subdir = os.path.dirname(_file).split("/")[1]
        filename = Path(_file).name
        os.makedirs(DATA_DIR / subdir, exist_ok=True)
        dest_path = (DATA_DIR / subdir / filename).resolve()
        if not dest_path.is_file():
            print("Downloading ", file_url, "->", str(dest_path))
            request.urlretrieve(file_url, dest_path)


if __name__ == "__main__":
    download_test_data()
