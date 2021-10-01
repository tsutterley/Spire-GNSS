"""
A tide prediction toolkit for Python
====================================

pyTMD contains Python tools for reading OTIS, GOT and FES formatted tidal
solutions to predict ocean and load tides

The package works using scientific Python packages (numpy, scipy and pyproj)
combined with data storage in netCDF4 and HDF5 and mapping using
matplotlib and cartopy

Documentation is available at https://pytmd.readthedocs.io
"""
import spire_toolkit.time
import spire_toolkit.utilities
import spire_toolkit.version

# get semantic version from setuptools-scm
__version__ = spire_toolkit.version.version
