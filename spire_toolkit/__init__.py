"""
A Spire GNSS toolkit for Python
===============================

spire_toolkit contains Python tools for obtaining and working with
elevation data from Spire GNSS grazing angle altimetry

The package works using Python packages (numpy, scipy, pyproj)
combined with data storage in netCDF4, and mapping with
matplotlib and cartopy

It aims to be a simple and efficient solution for using data from
Spire GNSS grazing angle altimetry and to support its science
applications

Documentation is available at https://spire-gnss.readthedocs.io
"""
import spire_toolkit.time
import spire_toolkit.utilities
import spire_toolkit.version

# get semantic version from setuptools-scm
__version__ = spire_toolkit.version.version
