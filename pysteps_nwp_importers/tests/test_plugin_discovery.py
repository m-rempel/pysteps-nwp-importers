import pytest
from pysteps.io import interface


new_importers = ["import_bom_nwp", "import_knmi_nwp", "import_rmi_nwp"]


@pytest.mark.parametrize("importer_name", new_importers)
def test_importers_discovery(importer_name):
    """Check if the importer plugin is correctly detected by pysteps. For this,
    the tests should be ran on the installed version of the plugin (and not
    against the plugin sources).
    """
    assert importer_name.replace("import_", "") in interface._importer_methods
