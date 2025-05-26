pysteps-nwp-importers
=====================

Pysteps plugin to import a variety of NWP rainfall forecasts. This plugin was created with Cookiecutter and makes the NWP importers findable for the Pysteps importer functionality after installation.


License
=======
* BSD license


Documentation
=============

NWP rainfall forecast importer using pysteps utilities. This plugin currently provides importers for:

* KNMI HARMONIE NWP forecasts
* BoM NWP forecasts
* RMI NWP forecasts
* DWD NWP forecasts


Installation
============

The latest development version of pysteps_nwp_importers can be installed using
pip by running in a terminal::

    pip install git+https://github.com/pySTEPS/pysteps-nwp-importers

Test the plugin
===============

This plugin comes with a tester, which can also be used to test whether the plugin is correctly hooked up to pysteps.

Install pytest and run the tests with::

	pip install pytest
	pytest -v --tb=line

Credits
=======

This package was created with Cookiecutter_ and the `cookiecutter-pysteps-plugin`_ project template. 
The `cookiecutter-pysteps-plugin`_ template was adapted from the cookiecutter-pypackage_
template.

.. _cookiecutter-pypackage: https://github.com/audreyfeldroy/cookiecutter-pypackage
.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`cookiecutter-pysteps-plugin`: https://github.com/pysteps/cookiecutter-pysteps-plugin