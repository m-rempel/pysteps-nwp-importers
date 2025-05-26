import os
from pathlib import Path
from urllib import request
from urllib.parse import urljoin

DATA_DIR = Path(__file__).parent / "../tests/data"


def download_test_data():
    root_url = "https://github.com/pySTEPS/pysteps-data/raw/master/nwp/"

    files_to_download = (
        "bom/2020/10/31/20201031_0000_regrid_short.nc",
        "knmi/2018/09/05/20180905_0600_Pforecast_Harmonie.nc",
        "rmi/2021/07/04/ao13_2021070412_native_5min.nc",
    )

    for _file in files_to_download:
        file_url = urljoin(root_url, _file)
        subdir = os.path.dirname(_file).split("/")[0]
        filename = Path(_file).name
        os.makedirs(DATA_DIR / subdir, exist_ok=True)
        dest_path = (DATA_DIR / subdir / filename).resolve()
        if not dest_path.is_file():
            print("Downloading ", file_url, "->", str(dest_path))
            request.urlretrieve(file_url, dest_path)


if __name__ == "__main__":
    download_test_data()
