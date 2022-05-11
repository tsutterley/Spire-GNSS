#!/usr/bin/env python
u"""
combine_spire_gnss_tiles.py
Written by Tyler Sutterley (05/2022)
Combines tile files of Spire GNSS surface fit data

COMMAND LINE OPTIONS:
    --help: list the command line options
    -D X, --directory X: directory with spire data
    -O X, --output X: Name and path of output file
    -H X, --hemisphere X: Region of interest to run
    -S X, --spacing X: Input tile grid spacing
    -s X, --subset X: Subset grid spacing for output
    -B X, --bounds X: valid range of tiles to read [xmin,xmax,ymin,ymax]
    -d, --dem: Input Digital Elevation Model file
    -o, --order: Order of surface polynomial fit
    -V, --verbose: Verbose output of run
    -M X, --mode X: Permissions mode of the directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    h5py: Python interface for Hierarchal Data Format 5 (HDF5)
        https://www.h5py.org/
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/

PROGRAM DEPENDENCIES:
    spatial.py: utilities for reading and writing spatial data
    time.py: Utilities for calculating time operations

UPDATE HISTORY:
    Updated 05/2022: use argparse descriptions within sphinx documentation
    Written 04/2022
"""
import sys
import os
import re
import h5py
import time
import pyproj
import netCDF4
import logging
import argparse
import warnings
import numpy as np
import pyTMD.time
import pyTMD.utilities
# from pyTMD.spatial import scale_areas
warnings.filterwarnings("ignore")

#-- PURPOSE: attempt to open an HDF5 file and wait if already open
def multiprocess_h5py(filename, *args, **kwargs):
    while True:
        try:
            fileID = h5py.File(filename, *args, **kwargs)
            break
        except (IOError, OSError, PermissionError) as e:
            time.sleep(0.2)
    return fileID

#-- PURPOSE: convert time from delta seconds into Julian and year-decimal
def convert_delta_time(delta_time, gps_epoch=0.0):
    """
    converts GPS delta_times into into Julian and year-decimal

    Arguments
    ---------
    delta_time: seconds since gps_epoch

    Keyword arguments
    -----------------
    gps_epoch: seconds between delta_time and GPS epoch (1980-01-06T00:00:00)

    Returns
    -------
    julian: time in Julian days
    decimal: time in year-decimal
    """
    #-- convert to array if single value
    delta_time = np.atleast_1d(delta_time)
    #-- calculate gps time from delta_time
    gps_seconds = gps_epoch + delta_time
    time_leaps = pyTMD.time.count_leap_seconds(gps_seconds)
    #-- calculate Julian time (UTC) by converting to MJD and then adding offset
    time_julian = 2400000.5 + pyTMD.time.convert_delta_time(
        gps_seconds - time_leaps, epoch1=(1980,1,6,0,0,0),
        epoch2=(1858,11,17,0,0,0), scale=1.0/86400.0)
    #-- convert to calendar date
    Y,M,D,h,m,s = pyTMD.time.convert_julian(time_julian,FORMAT='tuple')
    #-- calculate year-decimal time (UTC)
    time_decimal = pyTMD.time.convert_calendar_decimal(Y,M,day=D,
        hour=h,minute=m,second=s)
    #-- return both the Julian and year-decimal formatted dates
    return dict(julian=np.squeeze(time_julian),decimal=np.squeeze(time_decimal))

def scale_areas(lat, flat=1.0/298.257223563, ref=70.0):
    """
    Calculates area scaling factors for a polar stereographic projection
    including special case of at the exact pole

    Inputs:
        lat: latitude (degrees north)

    Options:
        flat: ellipsoidal flattening (default: WGS84)
        ref: reference latitude (true scale latitude)

    Returns:
        scale: area scaling factors at input latitudes

    References:
        Snyder, J P (1982) Map Projections used by the U.S. Geological Survey
            Forward formulas for the ellipsoid.  Geological Survey Bulletin
            1532, U.S. Government Printing Office.
        JPL Technical Memorandum 3349-85-101
    """
    #-- convert latitude from degrees to positive radians
    theta = np.abs(lat)*np.pi/180.0
    #-- convert reference latitude from degrees to positive radians
    theta_ref = np.abs(ref)*np.pi/180.0
    #-- square of the eccentricity of the ellipsoid
    #-- ecc2 = (1-b**2/a**2) = 2.0*flat - flat^2
    ecc2 = 2.0*flat - flat**2
    #-- eccentricity of the ellipsoid
    ecc = np.sqrt(ecc2)
    #-- calculate ratio at input latitudes
    m = np.cos(theta)/np.sqrt(1.0 - ecc2*np.sin(theta)**2)
    t = np.tan(np.pi/4.0 - theta/2.0)/((1.0 - ecc*np.sin(theta)) / \
        (1.0 + ecc*np.sin(theta)))**(ecc/2.0)
    #-- calculate ratio at reference latitude
    mref = np.cos(theta_ref)/np.sqrt(1.0 - ecc2*np.sin(theta_ref)**2)
    tref = np.tan(np.pi/4.0 - theta_ref/2.0)/((1.0 - ecc*np.sin(theta_ref)) / \
        (1.0 + ecc*np.sin(theta_ref)))**(ecc/2.0)
    #-- distance scaling
    k = (mref/m)*(t/tref)
    kp = 0.5*mref*np.sqrt(((1.0+ecc)**(1.0+ecc))*((1.0-ecc)**(1.0-ecc)))/tref
    #-- area scaling
    scale = np.where(np.isclose(theta,np.pi/2.0),1.0/(kp**2),1.0/(k**2))
    return scale

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Combine tiles of Spire GNSS surface fit data
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    #-- command line parameters
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Working data directory')
    #-- output combined file
    parser.add_argument('--output','-O',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Name and path of output file')
    #-- region of interest to run
    parser.add_argument('--hemisphere','-H',
        type=str, default='S', choices=('N','S'),
        help='Region of interest to run')
    #-- output grid spacing
    parser.add_argument('--spacing','-S',
        type=float, default=80e3, nargs='?',
        help='Output grid spacing')
    #-- subset grid spacing
    parser.add_argument('--subset','-s',
        type=float, default=1e3, nargs='?',
        help='Subset grid spacing')
    #-- bounds of output mosaic
    parser.add_argument('--bounds','-B', type=float,
        nargs=4, default=[-np.inf,np.inf,-np.inf,np.inf],
        metavar=('xmin','xmax','ymin','ymax'),
        help='valid range of tiles to read')
    #-- input digital elevation model for filtering and fit
    parser.add_argument('--dem','-d',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Input Digital Elevation Model file')
    #-- order of surface polynomial fit
    parser.add_argument('--order','-o',
        type=int, default=0,
        help='Order of surface polynomial fit')
    #-- verbose will output information about each output file
    parser.add_argument('--verbose','-V',
        action='count', default=0,
        help='Verbose output of processing run')
    #-- permissions mode of the directories and files (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files')
    #-- return the parser
    return parser

#-- PURPOSE: combine tiles of Spire GNSS surface fit data
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    #-- index directory for hemisphere
    index_directory = 'north' if (args.hemisphere == 'N') else 'south'
    DIRECTORY = os.path.join(args.directory, index_directory, 'fit')
    logging.info('Tile Fit Directory: {0}'.format(DIRECTORY))

    #-- regular expression pattern for tile files
    R1 = re.compile(r'E([-+]?\d+)_N([-+]?\d+)', re.VERBOSE)
    #-- output file format
    output_format = '{0}_{1}_{2}_{3:0.0f}km.h5'
    #-- invalid values for floating point variables
    fill_value = -9999.0

    #-- pyproj transformer for converting to polar stereographic
    EPSG = dict(N=3413, S=3031)[args.hemisphere]
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- dictionary of coordinate reference system variables
    cs_to_cf = crs2.cs_to_cf()
    crs_to_dict = crs2.to_dict()
    #-- flattening and standard parallel of datum and projection
    crs_to_cf = crs2.to_cf()

    #-- valid range of tiles
    xmin,xmax,ymin,ymax = args.bounds
    #-- dimensions of all tiles
    dimensions = [None,None]
    #-- find tile files within index directory
    all_tile_files = [f for f in os.listdir(DIRECTORY) if R1.match(f)]
    #-- reduce tiles to range
    tile_files = []
    for f in all_tile_files:
        xc,yc = [float(item)*1.e3 for item in R1.search(f).groups()]
        if ((xc >= xmin) and (xc <= xmax) & (yc >= ymin) and (yc <= ymax)):
            tile_files.append(f)
    #-- print number of valid tiles
    n_files = len(tile_files)
    logging.info('Number of tiles: {0:d}'.format(n_files))

    #-- broadcast spacing to dimensions if using uniform
    SPACING = np.broadcast_to(np.atleast_1d(args.spacing),(2,))
    #-- subset grid spacing
    SUBSET = np.broadcast_to(np.atleast_1d(args.subset),(2,))
    #-- subset grid dimensions
    nx = int(SPACING[0]//SUBSET[0])
    ny = int(SPACING[1]//SUBSET[1])
    area_km2 = SUBSET[0]*SUBSET[1]/1e6
    #-- calculate y dimensions with available extents
    dimensions[0] = np.int64((ymax - ymin)/SUBSET[1]) + ny
    #-- calculate x dimensions with available extents
    dimensions[1] = np.int64((xmax - xmin)/SUBSET[0]) + nx
    #-- print dimensions of grid
    logging.info('Grid Spacing: {0:0.0f} {1:0.0f}'.format(*SPACING))
    logging.info('Grid Dimensions: {0:d} {1:d}'.format(*dimensions))
    #-- print dimensions of subset grid
    logging.info('Subset Spacing: {0:0.0f} {1:0.0f}'.format(*SUBSET))
    logging.info('Subset Dimensions: {0:d} {1:d}'.format(nx,ny))
    #-- histogram parameters
    w = 0.05
    b = np.arange(0,8+w,w)
    nbins = len(b) - 1

    #-- output variable attributes
    attributes = {}
    #-- projection attributes
    attributes['Polar_Stereographic'] = {}
    #-- add projection attributes
    attributes['Polar_Stereographic']['standard_name'] = 'Polar_Stereographic'
    attributes['Polar_Stereographic']['spatial_epsg'] = crs2.to_epsg()
    attributes['Polar_Stereographic']['spatial_ref'] = crs2.to_wkt()
    attributes['Polar_Stereographic']['proj4_params'] = crs2.to_proj4()
    attributes['Polar_Stereographic']['latitude_of_projection_origin'] = \
        crs_to_dict['lat_0']
    for att_name,att_val in crs_to_cf.items():
        attributes['Polar_Stereographic'][att_name] = att_val
    #-- x and y
    attributes['x'],attributes['y'] = ({},{})
    for att_name in ['long_name','standard_name','units']:
        attributes['x'][att_name] = cs_to_cf[0][att_name]
        attributes['y'][att_name] = cs_to_cf[1][att_name]
    #-- histogram bins
    attributes['bin'] = {}
    attributes['bin']['units'] = 'm'
    attributes['bin']['long_name'] = 'Bin'
    attributes['bin']['description'] = ('Histogram bins of elevation '
        'over {0:0.0f}km tiles').format(area_km2)
    #-- elevation
    attributes['h_mean'] = {}
    attributes['h_mean']['units'] = 'm'
    attributes['h_mean']['long_name'] = 'Elevation'
    attributes['h_mean']['description'] = ('Average elevation '
        'estimated over {0:0.0f}km tiles').format(area_km2)
    attributes['h_mean']['grid_mapping'] = 'Polar_Stereographic'
    attributes['h_mean']['poly_order'] = args.order
    attributes['h_mean']['coordinates'] = 'y x'
    attributes['h_mean']['_FillValue'] = fill_value
    #-- elevation sigma
    attributes['h_sigma'] = {}
    attributes['h_sigma']['units'] = 'm'
    attributes['h_sigma']['long_name'] = 'Elevation uncertainty'
    attributes['h_sigma']['description'] = ('RMS of elevation '
        'uncertainty estimated over {0:0.0f}km tiles').format(area_km2)
    attributes['h_sigma']['grid_mapping'] = 'Polar_Stereographic'
    attributes['h_mean']['poly_order'] = args.order
    attributes['h_sigma']['coordinates'] = 'y x'
    attributes['h_sigma']['_FillValue'] = fill_value
    #-- elevation
    attributes['dem_h'] = {}
    attributes['dem_h']['units'] = 'm'
    attributes['dem_h']['long_name'] = 'Model Elevation'
    attributes['dem_h']['description'] = ('Average elevation '
        'estimated over {0:0.0f}km tiles').format(area_km2)
    attributes['dem_h']['grid_mapping'] = 'Polar_Stereographic'
    attributes['dem_h']['source'] = os.path.basename(args.dem)
    attributes['dem_h']['coordinates'] = 'y x'
    attributes['dem_h']['_FillValue'] = fill_value
    #-- grid cell areas
    attributes['cell_area'] = {}
    attributes['cell_area']['units'] = 'km^2'
    attributes['cell_area']['long_name'] = 'cell_area'
    attributes['cell_area']['description'] = ('Grid cell areas accounting for '
        'polar stereographic distortion')
    attributes['cell_area']['reference'] = 'Snyder, J P (1982)'
    attributes['cell_area']['grid_mapping'] = 'Polar_Stereographic'
    attributes['cell_area']['coordinates'] = 'y x'
    attributes['cell_area']['_FillValue'] = fill_value
    #-- point count
    attributes['data_count'] = {}
    attributes['data_count']['units'] = '1'
    attributes['data_count']['long_name'] = 'Count'
    attributes['data_count']['description'] = ('Point count within '
        '{0:0.0f}km tiles').format(area_km2)
    attributes['data_count']['grid_mapping'] = 'Polar_Stereographic'
    attributes['data_count']['coordinates'] = 'y x'
    attributes['data_count']['_FillValue'] = 0
    #-- error histogram
    attributes['histogram'] = {}
    attributes['histogram']['units'] = '1'
    attributes['histogram']['long_name'] = 'Histogram'
    attributes['histogram']['description'] = ('Histograms of elevation '
        'residuals within {0:0.0f}km tiles').format(area_km2)
    attributes['histogram']['grid_mapping'] = 'Polar_Stereographic'
    attributes['histogram']['coordinates'] = 'y x bin'
    attributes['histogram']['_FillValue'] = 0

    #-- dictionary with output variables
    output = {}
    #-- floating point variables
    variables = []
    variables.append('h_mean')
    variables.append('h_sigma')
    variables.append('dem_h')
    variables.append('cell_area')
    for key in variables:
        output[key] = np.ma.zeros((dimensions[0],dimensions[1]),
            fill_value=fill_value, dtype=np.float32)
        output[key].mask = np.ones((dimensions[0],dimensions[1]),
            dtype=bool)
    #-- number of spire GNSS points
    output['data_count'] = np.zeros((dimensions[0],dimensions[1]),
        dtype=np.uint32)
    #-- combined elevation error histograms
    output['histogram'] = np.zeros((dimensions[0],dimensions[1],nbins),
        dtype=np.uint32)
    #-- centers of each histogram bin
    output['bin'] = (b[1:] + b[0:-1])/2.0
    #-- projection variable
    output['Polar_Stereographic'] = np.empty((),dtype=np.byte)
    #-- grid coordinates
    output['x'] = np.linspace(xmin,xmax,dimensions[1])
    output['y'] = np.linspace(ymin,ymax,dimensions[0])

    #-- for each valid tile
    for iteration in range(n_files):
        #-- open the HDF5 tile file
        FILE1 = os.path.join(DIRECTORY,tile_files[iteration])
        f1 = h5py.File(FILE1, 'r')
        #-- get file-level attributes
        xc = f1.attrs['x_center']
        yc = f1.attrs['y_center']
        logging.debug('Tile File: {0}'.format(tile_files[iteration]))
        logging.debug('Tile Center: {0:0.0f} {1:0.0f}'.format(xc,yc))
        #-- subset grid coordinates
        xsub = f1['x'][:]
        ysub = f1['y'][:]
        #-- image coordinates of subset grid
        IMx = np.array((xsub-xmin)/SUBSET[0],dtype=np.int64)
        IMy = np.array((ysub-ymin)/SUBSET[1],dtype=np.int64)
        ix,iy = slice(IMx[0],IMx[-1]), slice(IMy[0],IMy[-1])
        try:
            #-- add image to output mosaic
            for key in variables:
                output[key].data[iy,ix] = f1[key][:,:].copy()
            #-- add data count to output mosaic
            output['data_count'][iy,ix] += f1['data_count'][:,:].copy()
            #-- add histogram to output mosaic
            output['histogram'][iy,ix,:] = f1['histogram'][...].copy()
        except:
            pass
        else:
            #-- update mosaic mask
            for key in variables:
                fv = f1[key].fillvalue
                output[key].mask[iy,ix] = (f1[key][:,:] == fv)
        #-- close the tile file
        f1.close()

    #-- replace invalid points with fill value
    for key in variables:
        output[key].data[output[key].mask] = output[key].fill_value

    #-- output HDF5 filename
    if args.output:
        #-- use defined output file name and path
        FILE2 = os.path.expanduser(args.output)
    else:
        #-- use default output file name and path
        fargs = ('spire_gnss', 'L2_grzAlt', 'h_mean', area_km2)
        FILE2 = os.path.join(DIRECTORY,output_format.format(*fargs))
    #-- recursively create output directory if not existent
    if not os.access(os.path.dirname(FILE2),os.F_OK):
        os.makedirs(os.path.dirname(FILE2), mode=args.mode)
    #-- opening HDF5 file for writing
    fileID = h5py.File(FILE2, 'w')
    #-- Defining the HDF5 dataset variables
    h5 = {}
    for key,val in output.items():
        if '_FillValue' in attributes[key].keys():
            h5[key] = fileID.create_dataset(key, val.shape, data=val,
                dtype=val.dtype, fillvalue=attributes[key]['_FillValue'],
                compression='gzip')
        elif val.shape:
            h5[key] = fileID.create_dataset(key, val.shape, data=val,
                dtype=val.dtype, compression='gzip')
        else:
            h5[key] = fileID.create_dataset(key, val.shape,
                dtype=val.dtype)
        #-- Defining attributes for variable
        for att_name,att_val in attributes[key].items():
            h5[key].attrs[att_name] = att_val
    #-- add attribute for date created
    today = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    fileID.attrs['date_created'] = today
    #-- add file-level attributes
    fileID.attrs['GDAL_AREA_OR_POINT'] = 'Area'
    fileID.attrs['source'] = 'Spacecraft'
    #-- convert regional coordinates to latitude and longitude
    xgrid,ygrid = np.meshgrid(output['x'],output['y'])
    direction = pyproj.enums.TransformDirection.INVERSE
    lon,lat = transformer.transform(xgrid,ygrid,direction=direction)
    ltmn,ltmx = np.min(lat),np.max(lat)
    lnmn,lnmx = np.min(lon),np.max(lon)
    #-- set attribute for total valid data count
    fileID.attrs['valid_count'] = np.sum(output['data_count'])
    fileID.attrs['geospatial_h_min'] = np.min(output['h_mean'])
    fileID.attrs['geospatial_h_max'] = np.max(output['h_mean'])
    fileID.attrs['geospatial_h_units'] = 'metres'
    #-- add geospatial attributes
    fileID.attrs['geospatial_lat_min'] = ltmn
    fileID.attrs['geospatial_lat_max'] = ltmx
    fileID.attrs['geospatial_lon_min'] = lnmn
    fileID.attrs['geospatial_lon_max'] = lnmx
    fileID.attrs['geospatial_lat_units'] = "degrees_north"
    fileID.attrs['geospatial_lon_units'] = "degrees_east"
    fileID.attrs['geospatial_ellipsoid'] = "WGS84"
    #-- Output HDF5 structure information
    logging.info(FILE2)
    logging.info(list(fileID.keys()))
    #-- Closing the HDF5 file
    fileID.close()
    #-- change the permissions mode
    os.chmod(FILE2, args.mode)

#-- run main program
if __name__ == '__main__':
    main()
