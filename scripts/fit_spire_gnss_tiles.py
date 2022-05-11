#!/usr/bin/env python
u"""
fit_spire_gnss_tiles.py
Written by Tyler Sutterley (05/2022)
Fits elevation surfaces to tile files of Spire GNSS data

INPUTS:
    infile: Spire GNSS tile file

COMMAND LINE OPTIONS:
    --help: list the command line options
    -O X, --output X: Name and path of output file
    -H X, --hemisphere X: Region of interest to run
    -S X, --spacing X: Input tile grid spacing
    -s X, --subset X: Subset grid spacing for output
    -f, --filter: Filter elevations using median statistics
    -d, --dem: Input Digital Elevation Model file
    -o, --order: Order of surface polynomial fit
    -V, --verbose: Verbose output of run
    -M X, --mode X: Permissions mode of the directories and files

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/
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
    Updated 11/2021: output merged tile file with filenames
        adjust tiling to index by center coordinates
        wait if merged netCDF4 tile file is unavailable
    Written 10/2021
"""
import sys
import os
import re
import h5py
import time
import pyproj
import logging
import netCDF4
import argparse
import warnings
import numpy as np
import pyTMD.time
import pyTMD.spatial
import pyTMD.utilities
warnings.filterwarnings("ignore")

#-- PURPOSE: attempt to open an netCDF4 file and wait if already open
def multiprocess_netCDF4(filename, *args, **kwargs):
    while True:
        try:
            fileID = netCDF4.Dataset(filename, *args, **kwargs)
            break
        except (IOError, OSError, PermissionError) as e:
            time.sleep(1)
    return fileID

#-- PURPOSE: convert time from delta seconds into Julian and year-decimal
def convert_delta_time(delta_time, gps_epoch=0.0):
    """
    converts spire delta_times into into Julian and year-decimal

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

#-- PURPOSE: check if h is valid by checking the interquartile range from
#-- Pritchard (2009) and the robust dispersion estimator (RDE) from Smith (2017)
def dh_filter(h0):
    #-- validate that heights are valid
    isfinite = np.nonzero(np.isfinite(h0))
    IQR_valid = np.zeros_like(h0,dtype=bool)
    RDE_valid = np.zeros_like(h0,dtype=bool)
    #-- check if all points are invalid
    if not np.count_nonzero(np.isfinite(h0)):
        return (IQR_valid,RDE_valid)
    #-- calculate percentiles for IQR, MDE and median
    #-- IQR: first and third quartiles (25th and 75th percentiles)
    #-- MDE: 16th and 84th percentiles
    #-- median: 50th percentile
    Q1,Q3,P16,P84,MEDIAN = np.percentile(h0[isfinite],[25,75,16,84,50])
    #-- calculate interquartile range
    IQR = Q3 - Q1
    #-- calculate robust dispersion estimator (RDE)
    RDE = P84 - P16
    #-- IQR pass: dh/dt of point-(median value) is within 75% of IQR
    IQR_valid[isfinite] = (np.abs(h0[isfinite]-MEDIAN) <= (0.75*IQR))
    #-- RDE pass: dh/dt of point-(median value) is within 50% of P84-P16
    RDE_valid[isfinite] = (np.abs(h0[isfinite]-MEDIAN) <= (0.50*RDE))
    #-- return the valid flags
    return (IQR_valid,RDE_valid)

#-- PURPOSE: fit a surface polynomial to height measurements
def fit_poly_surface(x, y, h, order=1, centroid=None):
    #-- remove singleton dimensions from input variables
    x = np.squeeze(x)
    y = np.squeeze(y)
    h = np.squeeze(h)
    nmax = len(h)
    #-- calculate x and y relative to centroid point
    rel_x = x - centroid['x']
    rel_y = y - centroid['y']
    #-- Constant Term
    P_x0 = np.ones((nmax))
    #-- Surface design matrix
    if (order == 0):
        DMAT = np.transpose([P_x0])
    elif (order == 1):#-- planar surface fit
        DMAT = np.transpose([P_x0, rel_x, rel_y])
    elif (order == 2):#-- quadratic surface fit
        DMAT = np.transpose([P_x0, rel_x, rel_y, rel_x**2, rel_y**2, rel_x*rel_y])
    elif (order == 3):#-- cubic surface fit
        DMAT = np.transpose([P_x0, rel_x, rel_y, rel_x**2, rel_y**2, rel_x*rel_y,
            rel_x**3, rel_y**3, (rel_x)*(rel_y**2), (rel_x**2)*(rel_y)])
    elif (order == 4):#-- quartic surface fit
        DMAT = np.transpose([P_x0, rel_x, rel_y, rel_x**2, rel_y**2, rel_x*rel_y,
            rel_x**3, rel_y**3, (rel_x)*(rel_y**2), (rel_x**2)*(rel_y),
            rel_x**4, rel_y**4, (rel_x)*(rel_y**3), (rel_x**3)*(rel_y),
            (rel_x**2)*(rel_y**2)])
    #-- Standard Least-Squares fitting (the [0] denotes coefficients output)
    beta_mat = np.linalg.lstsq(DMAT,h,rcond=-1)[0]
    n_terms = len(beta_mat)
    #-- modelled surface with all parameters
    mod = np.dot(DMAT,beta_mat)
    #-- residual of fit
    res = h - mod
    #-- nu = Degrees of Freedom
    nu = nmax - n_terms
    #-- Mean square error
    #-- MSE = (1/nu)*sum((Y-X*B)**2)
    MSE = np.dot(np.transpose(h - np.dot(DMAT,beta_mat)),
        (h - np.dot(DMAT,beta_mat)))/nu
    #-- return the fit parameters
    return {'beta':beta_mat, 'model':mod, 'residual': res, 'MSE':MSE, 'DOF':nu}

#-- PURPOSE: create argument parser
def arguments():
    parser = argparse.ArgumentParser(
        description="""Fits elevation surfaces to tile files of Spire GNSS data
            """,
        fromfile_prefix_chars="@"
    )
    parser.convert_arg_line_to_args = pyTMD.utilities.convert_arg_line_to_args
    #-- command line parameters
    parser.add_argument('infile',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        help='Input Spire GNSS tile file')
    #-- output combined spire tile file
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
        help='Input tile grid spacing')
    #-- subset grid spacing
    parser.add_argument('--subset','-s',
        type=float, default=1e3, nargs='?',
        help='Subset grid spacing for output')
    #-- filter elevations using median statistics
    parser.add_argument('--filter','-f',
        default=False, action='store_true',
        help='Filter elevations using median statistics')
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

#-- PURPOSE: fit elevation surfaces to tiles of Spire GNSS data
def main():
    #-- Read the system arguments listed after the program
    parser = arguments()
    args,_ = parser.parse_known_args()

    #-- create logger
    loglevels = [logging.CRITICAL,logging.INFO,logging.DEBUG]
    logging.basicConfig(level=loglevels[args.verbose])

    #-- index directory for hemisphere
    DIRECTORY = os.path.dirname(args.infile)
    file_directory = os.path.dirname(DIRECTORY)
    logging.info('Tile directory: {0}'.format(DIRECTORY))

    #-- regular expression pattern for tile files
    R1 = re.compile(r'E([-+]?\d+)_N([-+]?\d+)', re.VERBOSE)
    #-- regular expression pattern for spire GNSS files
    R2 = re.compile(r'(spire_gnss-r)_(L\d+)_(.*?)_(v\d+\.\d+)_'
        r'(\d{4})-(\d{2})-(\d{2})T(\d{2})-(\d{2})-(\d{2})_(.*?)\.nc$')
    #-- invalid values for floating point variables
    fill_value = -9999.0
    #-- count threshold for valid fits
    count_threshold = 5

    #-- pyproj transformer for converting to polar stereographic
    EPSG = dict(N=3413,S=3031)[args.hemisphere]
    crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(EPSG))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    #-- dictionary of coordinate reference system variables
    cs_to_cf = crs2.cs_to_cf()
    crs_to_dict = crs2.to_dict()
    #-- flattening and standard parallel of datum and projection
    crs_to_cf = crs2.to_cf()
    flattening = 1.0/crs_to_cf['inverse_flattening']
    standard_parallel = crs_to_cf['standard_parallel']

    #-- broadcast spacing to dimensions if using uniform
    SPACING = np.broadcast_to(np.atleast_1d(args.spacing),(2,))
    #-- print dimensions of grid
    logging.info('Grid Spacing: {0:0.0f} {1:0.0f}'.format(*SPACING))

    #-- subset grid spacing
    SUBSET = np.broadcast_to(np.atleast_1d(args.subset),(2,))
    area_km2 = SUBSET[0]*SUBSET[1]/1e6
    nx = int(SPACING[0]//SUBSET[0]) + 1
    ny = int(SPACING[1]//SUBSET[1]) + 1
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
        output[key] = np.ma.zeros((ny,nx), fill_value=fill_value,
            dtype=np.float32)
        output[key].mask = np.ones((ny,nx),dtype=bool)
    #-- number of spire GNSS points
    output['data_count'] = np.zeros((ny,nx),dtype=np.uint32)
    #-- combined elevation error histograms
    output['histogram'] = np.zeros((ny,nx,nbins),dtype=np.uint32)
    #-- centers of each histogram bin
    output['bin'] = (b[1:] + b[0:-1])/2.0
    #-- projection variable
    output['Polar_Stereographic'] = np.empty((),dtype=np.byte)

    #-- read DEM subsetted to tile
    if args.dem is not None:
        DEM = pyTMD.spatial.from_geotiff(args.dem)
        dx, dy = DEM['attributes']['spacing']
        xmin, xmax, ymin, ymax = DEM['attributes']['extent']
        fv = DEM['attributes']['data']['_FillValue']

    #-- open the HDF5 tile file
    f1 = netCDF4.Dataset(args.infile, 'r')
    #-- get file-level attributes
    xc = f1.x_center
    yc = f1.y_center
    logging.debug('Tile Center: {0:0.0f} {1:0.0f}'.format(xc,yc))
    #-- x and y coordinate for subset image
    xsub = np.linspace(xc-SPACING[0]/2, xc+SPACING[0]/2, nx)
    ysub = np.linspace(yc-SPACING[1]/2, yc+SPACING[1]/2, ny)
    #-- calculate x and y arrays of output grid
    output['x'] = xsub + SUBSET[0]/2.0
    output['y'] = ysub + SUBSET[1]/2.0
    #-- meshgrid of coordinate centers
    gridx,gridy = np.meshgrid(output['x'], output['y'])
    #-- convert grid coordinates from polar stereographic
    direction = pyproj.enums.TransformDirection.INVERSE
    _,gridlat = transformer.transform(gridx,gridy,direction=direction)
    #-- get area weights for each latitude
    output['cell_area'][:,:] = pyTMD.spatial.scale_areas(gridlat,
        flat=flattening, ref=standard_parallel)

    #-- find Spire GNSS files within tile file (groups)
    GNSS_files = [f for f in f1.groups if R2.match(f)]
    n_files = len(GNSS_files)
    logging.debug('Number of files: {0:d}'.format(n_files))

    #-- calculate total number of points in tile
    N = 0
    for iteration in range(n_files):
        #-- spire GNSS file for iteration
        GNSS = GNSS_files[iteration]
        indices = f1[GNSS].variables['index'][:].copy()
        N += len(indices)

    #-- dictionary with data
    data = {}
    for k in ['height','spec_lon','spec_lat','spec_elevation','spec_azimuth',
        'time','ocean_tide','slant_bias']:
        data[k] = np.zeros((N))
    #-- counter variable for filling arrays
    c = 0
    #-- read each Spire GNSS file and estimate heights
    for iteration in range(n_files):
        #-- spire GNSS file for iteration
        GNSS = GNSS_files[iteration]
        FILE2 = os.path.join(file_directory,GNSS)
        #-- extract parameters from netCDF4 file
        MS,LV,PRD,VERS,YY,MM,DD,HH,MN,SS,AUX = R2.findall(GNSS).pop()
        #-- open GNSS file
        f2 = multiprocess_netCDF4(FILE2, 'r')
        #-- reference points and indices within tile
        indices = f1[GNSS]['index'][:].copy()
        n_pts = len(indices)
        #-- extract GNSS variables for indices
        for k in ['height','spec_lon','spec_lat',
            'spec_elevation','spec_azimuth',
            'time','ocean_tide','slant_bias']:
            temp = f2.variables[k][:].copy()
            #-- reduce to indices
            data[k][c:c+n_pts] = temp[indices]
        #-- invalid value for sigmas
        invalid = f2.variables['height']._FillValue
        #-- add to counter
        c += n_pts
        #-- close the GNSS file
        f2.close()
    #-- close the GNSS tile file
    f1.close()

    #-- replace fill values with nan or zero
    data['height'][data['height'] == invalid] = np.nan
    data['ocean_tide'][data['ocean_tide'] == invalid] = 0.0
    #-- convert each pair to polar stereographic
    x, y = transformer.transform(data['spec_lon'],
        data['spec_lat'])
    #-- reduce heights to high quality estimates
    mask = np.isfinite(data['height'])
    #-- use DEM elevation as baseline for fit
    h_mod = np.zeros_like(data['height'])
    if args.dem is not None:
        #-- initially set all values to fill value
        h_mod[:] = fv
        #-- calculate image coordinates for DEM
        IMx = np.array((x - xmin)//dx, dtype='i')
        IMy = np.array((y - ymax)//dy, dtype='i')
        #-- limit to valid range of image coordinates
        i = np.nonzero((IMx >= 0) & (IMx <= DEM['data'].shape[1]) &
            (IMy >= 0) & (IMy <= DEM['data'].shape[1]))
        h_mod[i] = DEM['data'][IMy[i],IMx[i]]
        mask &= (h_mod != fv)
    #-- calculate relative elevation
    dh = data['height'] - h_mod

    #-- for each image coordinate
    for i,ys in enumerate(ysub):
        for j,xs in enumerate(xsub):
            #-- find where data is within subset pixel
            subset = ((x >= xs) & (x < (xs+SUBSET[0])) &
                (y >= ys) & (y < (ys+SUBSET[1])))
            good = np.ones_like(mask)
            #-- use median statistics to filter elevations
            if args.filter:
                ind, = np.nonzero(mask & subset)
                IQR,RDE = dh_filter(dh[ind])
                good[ind] &= RDE
            #-- valid heights with spatial subset
            count = np.count_nonzero(mask & subset & good)
            valid, = np.nonzero(mask & subset & good)
            #-- check that there are some valid points
            if (count < count_threshold):
                continue
            #-- fit relative to center of grid cell
            centroid = {}
            centroid['x'] = (xs + SUBSET[0]/2.0)
            centroid['y'] = (ys + SUBSET[1]/2.0)
            #-- fit a polynomial surface
            fit = fit_poly_surface(x[valid], y[valid], dh[valid],
                order=args.order, centroid=centroid)
            if args.dem is not None:
                #-- calculate image coordinates for DEM
                XS = slice(int((xs - xmin)//dx),int((xs + SUBSET[0] - xmin)//dx))
                YS = slice(int((ys + SUBSET[1] - ymax)//dy),int((ys - ymax)//dy))
                output['dem_h'][i,j] = np.mean(DEM['data'][YS,XS])
            else:
                output['dem_h'][i,j] = 0.0
            #-- average elevation and elevation uncertainty
            output['h_mean'][i,j] = fit['beta'][0] + output['dem_h'][i,j]
            output['h_sigma'][i,j] = np.sqrt(fit['MSE'])
            #-- add to histograms
            h,_ = np.histogram(np.abs(fit['residual']), bins=b)
            output['histogram'][i,j,:] += h.astype(np.uint32)
            #-- add to data count
            output['data_count'][i,j] = np.copy(count)
            #-- update mask for variables
            for key in variables:
                output[key].mask[i,j] = False

    #-- replace invalid points with fill value
    for key in variables:
        output[key].data[output[key].mask] = output[key].fill_value

    #-- output HDF5 filename
    if args.output:
        #-- use defined output file name and path
        FILE4 = os.path.expanduser(args.output)
    else:
        #-- use default output file name and path
        fileBasename,_ = os.path.splitext(os.path.basename(args.infile))
        FILE4 = os.path.join(DIRECTORY,'fit','{0}.h5'.format(fileBasename))
    #-- recursively create output directory if not existent
    if not os.access(os.path.dirname(FILE4),os.F_OK):
        os.makedirs(os.path.dirname(FILE4), mode=args.mode)
    #-- opening HDF5 file for writing
    fileID = h5py.File(FILE4, 'w')
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
    #-- set attribute for tile center
    fileID.attrs['x_center'] = xc
    fileID.attrs['y_center'] = yc
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
    logging.info(FILE4)
    logging.info(list(fileID.keys()))
    #-- Closing the HDF5 file
    fileID.close()
    #-- change the permissions mode
    os.chmod(FILE4, args.mode)

#-- run main program
if __name__ == '__main__':
    main()
