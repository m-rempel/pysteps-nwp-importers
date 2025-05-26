from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

requirements = [
    "numpy",
    "pygrib",
    "pyproj",
    "xarray",
    "dask",
    "netCDF4",
    "pysteps>=1.4,<2.0",
]

setup_requirements = ["pytest-runner"]

test_requirements = ["pytest>=3"]

setup(
    author="pySTEPS developers",
    python_requires=">=3.7",
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    description="Pysteps plugin to import a variety of NWP rainfall forecasts",
    install_requires=requirements,
    license="BSD license",
    long_description=readme,
    test_suite="tests",
    tests_require=test_requirements,
    include_package_data=True,
    keywords=["pysteps_nwp_importers", "pysteps", "plugin", "importer"],
    name="pysteps_nwp_importers",
    packages=find_packages(),
    setup_requires=setup_requirements,
    entry_points={
        "pysteps.plugins.importers": [
            "import_knmi_nwp=pysteps_nwp_importers.importer_knmi_nwp:import_knmi_nwp",
            "import_bom_nwp=pysteps_nwp_importers.importer_bom_nwp:import_bom_nwp",
            "import_rmi_nwp=pysteps_nwp_importers.importer_rmi_nwp:import_rmi_nwp",
            "import_dwd_nwp=pysteps_nwp_importers.importer_dwd_nwp:import_dwd_nwp",
        ]
    },
    version="0.1",
    zip_safe=False,
)
