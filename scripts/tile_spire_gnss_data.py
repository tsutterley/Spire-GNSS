#!/usr/bin/env python
u"""
tile_spire_gnss_data.py
Written by Tyler Sutterley (05/2022)
Creates tile index files of Spire GNSS data

INPUTS:
    input_file: Spire GNSS data file

COMMAND LINE OPTIONS:
    --help: list the command line options
    -H X, --hemisphere X: Region of interest to run
    -S X, --spacing X: Output grid spacing
    -V, --verbose: Verbose output of run
    -M X, --mode X: Permissions mode of the directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://www.numpy.org
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Updated 11/2021: output merged tile file with filenames
        adjust tiling to index by center coordinates
        wait if merged netCDF4 tile file is unavailable
    Written 10/2021
"""
import sys
import os
import re
import time
import pyproj
import logging
import netCDF4
import argparse
import numpy as np

#-- PURPOSE: attempt to open an netCDF4 file and wait if already open
def multiprocess_netCDF4(filename, *args, **kwargs):
    while True:
        try:
            fileID = netCDF4.Dataset(filename, *args, **kwargs)
            break
        except (IOError, OSError, PermissionError) as e:
            time.sleep(1)
    return fileID

#-- PURPOSE: create tile index files of Spire GNSS data
def tile_spire_gnss_data(input_file,
    SPACING=None,
    HEM=None,
    VERBOSE=False,
    MODE=0o775):

    #-- create logger
    loglevel = logging.INFO if VERBOSE else logging.CRITICAL
    logging.basicConfig(level=loglevel)

    #-- index directory for hemisphere
    index_directory = 'north' if (HEM == 'N') else 'south'
    #-- output directory and index file
    DIRECTORY = os.path.dirname(input_file)
    BASENAME = os.path.basename(input_file)
    output_file = os.path.join(DIRECTORY, index_directory, BASENAME)
    #-- regular expression pattern for extracting data from files
    rx = re.compile(r'(spire_gnss-r)_(L\d+)_(.*?)_(v\d+\.\d+)_'
        r'(\d{4})-(\d{2})-(\d{2})T(\d{2})-(\d{2})-(\d{2})_(.*?)\.nc$')
    #-- extract parameters from netCDF4 file
    MS,LV,PRD,VERS,YY,MM,DD,HH,MN,SS,AUX = rx.findall(input_file).pop()

    #-- pyproj transformer for converting to polar stereographic
    EPSG = dict(N=3413,S=3031)
    SIGN = dict(N=1.0,S=-1.0)
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG[HEM]))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- dictionary of coordinate reference system variables
    cs_to_cf = crs2.cs_to_cf()
    crs_to_dict = crs2.to_dict()

    #-- attributes for each output item
    attributes = dict(x={},y={},index={},time={})
    #-- x and y
    for att_name in ['long_name','standard_name','units']:
        attributes['x'][att_name] = cs_to_cf[0][att_name]
        attributes['y'][att_name] = cs_to_cf[1][att_name]
    #-- index
    attributes['index']['long_name'] = 'Index'
    attributes['index']['grid_mapping'] = 'Polar_Stereographic'
    attributes['index']['units'] = '1'
    attributes['index']['coordinates'] = 'x y'
    #-- time
    attributes['time']['long_name'] = 'time'
    attributes['time']['standard_name'] = 'time'
    attributes['time']['description'] = 'Time of GPS measurement (GPS seconds)'
    attributes['time']['units'] = 'seconds since 1980-01-06T00:00:00Z'
    attributes['time']['calendar'] = 'standard'

    #-- create index directory for hemisphere
    if not os.access(os.path.join(DIRECTORY,index_directory),os.F_OK):
        os.makedirs(os.path.join(DIRECTORY,index_directory),
            mode=MODE, exist_ok=True)

    #-- track file progress
    logging.info(input_file)
    #-- open spire file
    with netCDF4.Dataset(input_file,'r') as f1:
        #-- extract latitude, longitude and time variables
        spec_lat = f1.variables['spec_lat'][:].copy()
        spec_lon = f1.variables['spec_lon'][:].copy()
        gps_time = f1.variables['time'][:].copy()

    #-- check that data is in hemisphere
    if not np.any(np.sign(spec_lat) == SIGN[HEM]):
        return
    #-- indices of points in hemisphere
    valid, = np.nonzero(np.sign(spec_lat) == SIGN[HEM])
    #-- convert latitude and longitude to regional projection
    x,y = transformer.transform(spec_lon,spec_lat)
    #-- large-scale tiles
    xtile = (x-0.5*SPACING)//SPACING
    ytile = (y-0.5*SPACING)//SPACING

    #-- open output index file
    f2 = netCDF4.Dataset(output_file,'w')
    f2.setncattr('featureType','trajectory')
    f2.setncattr('GDAL_AREA_OR_POINT','Point')
    f2.setncattr('Conventions','CF-1.6')
    f2.setncattr('time_type','GPS')
    today = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    f2.setncattr('date_created', today)
    #-- create projection variable
    nc = f2.createVariable('Polar_Stereographic',np.byte,())
    #-- add projection attributes
    nc.standard_name = 'Polar_Stereographic'
    nc.spatial_epsg = crs2.to_epsg()
    nc.spatial_ref = crs2.to_wkt()
    nc.proj4_params = crs2.to_proj4()
    nc.latitude_of_projection_origin = crs_to_dict['lat_0']
    for att_name,att_val in crs2.to_cf().items():
        nc.setncattr(att_name,att_val)
    #-- for each valid tile pair
    for xp,yp in set(zip(xtile[valid],ytile[valid])):
        #-- center of each tile (adjust due to integer truncation)
        xc = (xp+1)*SPACING
        yc = (yp+1)*SPACING
        #-- create group
        tile_group = 'E{0:0.0f}_N{1:0.0f}'.format(xc/1e3,yc/1e3)
        g2 = f2.createGroup(tile_group)
        #-- add group attributes
        g2.setncattr('x_center',xc)
        g2.setncattr('y_center',yc)
        g2.setncattr('spacing',SPACING)

        #-- create merged tile file if not existing
        tile_file = os.path.join(DIRECTORY,index_directory,
            '{0}.nc'.format(tile_group))
        clobber = 'a' if os.access(tile_file,os.F_OK) else 'w'
        #-- open output merged tile file
        f3 = multiprocess_netCDF4(tile_file, clobber)
        g3 = f3.createGroup(BASENAME)
        #-- add file-level variables and attributes
        if (clobber == 'w'):
            #-- create projection variable
            nc = f3.createVariable('Polar_Stereographic',np.byte,())
            #-- add projection attributes
            nc.standard_name = 'Polar_Stereographic'
            nc.spatial_epsg = crs2.to_epsg()
            nc.spatial_ref = crs2.to_wkt()
            for att_name,att_val in crs2.to_cf().items():
                nc.setncattr(att_name,att_val)
            #-- add file attributes
            f3.setncattr('featureType','trajectory')
            f3.setncattr('x_center',xc)
            f3.setncattr('y_center',yc)
            f3.setncattr('spacing',SPACING)
            f3.setncattr('GDAL_AREA_OR_POINT','Point')
            f3.setncattr('Conventions','CF-1.6')
            f3.setncattr('time_type','GPS')
            f3.setncattr('date_created', today)

        #-- indices of points within tile
        indices, = np.nonzero((xtile == xp) & (ytile == yp))
        #-- output variables for index file
        output = dict(x=x[indices],y=y[indices],
            time=gps_time[indices],index=indices.copy())
        #-- for each group
        for g in [g2,g3]:
            #-- Defining the netCDF4 dimensions
            g.createDimension('time', len(output['time']))
            for key,val in output.items():
                #-- Defining the netCDF4 variables
                nc = g.createVariable(key, val.dtype, ('time',))
                #-- filling the netCDF4 variables
                nc[:] = val.copy()
                #-- add variable attributes
                for att_name,att_val in attributes[key].items():
                    setattr(nc,att_name,att_val)
        #-- close the merged tile file
        f3.close()

    #-- Output netCDF4 structure information
    logging.info(list(f2.groups.keys()) + list(f2.variables.keys()))
    #-- close the output file
    f2.close()
    #-- change the permissions mode of the output file
    os.chmod(output_file, mode=MODE)

#-- PURPOSE: create arguments parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Creates tile index files of Spire GNSS data
            """
    )
    #-- command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)), nargs='+',
        help='Input Spire GNSS file')
    #-- region of interest to run
    parser.add_argument('--hemisphere','-H',
        type=str, default='N', choices=('N','S'),
        help='Hemisphere')
    #-- output grid spacing
    parser.add_argument('--spacing','-S',
        type=float, default=10e3,
        help='Output grid spacing')
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        default=False, action='store_true',
        help='Verbose output of run')
    #-- permissions mode of the directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    #-- return the parser
    return parser

#-- This is the main part of the program that calls the individual functions
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- run program for each product
    for FILE in args.infile:
        tile_spire_gnss_data(FILE,
            SPACING=args.spacing,
            HEM=args.hemisphere,
            VERBOSE=args.verbose,
            MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
