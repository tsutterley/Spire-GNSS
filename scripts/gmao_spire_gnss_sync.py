#!/usr/bin/env python
u"""
gmao_spire_gnss_sync.py
Written by Tyler Sutterley (10/2021)

Syncs Spire GNSS grazing angle altimetry data from the NASA
    Global Modeling and Assimilation Office (GMAO)

CALLING SEQUENCE:
    python gmao_spire_gnss_sync.py --user=<username>
    where <username> is your NASA GMAO Extranet credentials

COMMAND LINE OPTIONS:
    --help: list the command line options
    -U X, --user X: username for NASA GMAO Extranet Login
    -W X, --password X: password for NASA GMAO Extranet Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: working data directory
    -p X, --product X: Spire data products to sync
    -Y X, --year X: Years of Spire data to sync
    -P X, --np X: Number of processes to use in file downloads
    -t X, --timeout X: Timeout in seconds for blocking operations
    -l, --log: output log of files downloaded
    -M X, --mode X: Local permissions mode of the directories and files synced

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
    dateutil: powerful extensions to datetime
        https://dateutil.readthedocs.io/en/stable/
    lxml: Pythonic XML and HTML processing library using libxml2/libxslt
        https://lxml.de/
        https://github.com/lxml/lxml
    future: Compatibility layer between Python 2 and Python 3
        https://python-future.org/

PROGRAM DEPENDENCIES:
    utilities.py: download and management utilities for syncing files

UPDATE HISTORY:
    Updated 10/2021: using python logging for handling verbose output
    Written 10/2021
"""
from __future__ import print_function

import sys
import os
import io
import re
import time
import netrc
import getpass
import logging
import tarfile
import builtins
import argparse
import posixpath
import traceback
import multiprocessing as mp
import spire_toolkit.utilities

#-- PURPOSE: Syncs Spire GNSS grazing angle altimetry data
def gmao_spire_gnss_sync(DIRECTORY,
    PRODUCT=[],
    YEAR=None,
    PROCESSES=0,
    TIMEOUT=None,
    LOG=False,
    MODE=0o775):

    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: GMAO_Spire_GNSS_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'GMAO_Spire_GNSS_sync_{0}.log'.format(today)
        logging.basicConfig(filename=os.path.join(DIRECTORY,LOGFILE),
            level=logging.INFO)
        logging.info('GMAO Spire GNSS Sync Log ({0})'.format(today))
    else:
        #-- standard output (terminal output)
        logging.basicConfig(level=logging.INFO)

    #-- remote host for Spire GNSS data
    HOST = 'https://gmao.gsfc.nasa.gov'
    #-- regular expression pattern for finding tar files
    regex_pattern = r'(spire_gnss-r)_(L\d+).({0})_({1})(\d{{2}}).tar'
    regex_products = r'|'.join(PRODUCT) if PRODUCT else r'.*?'
    regex_years = r'|'.join([r'{0:4d}'.format(y) for y in YEAR])
    R1 = re.compile(regex_pattern.format(regex_products,regex_years))
    #-- open connection with GMAO extranet server at remote directory
    PATH = [HOST,'extranet','collab','spire_team']
    files = spire_toolkit.utilities.gmao_list(PATH, timeout=TIMEOUT,
        adddirlink="grazing_angle_L2", pattern=R1, sort=True)

    #-- sync in series if PROCESSES = 0
    if (PROCESSES == 0):
        #-- get each tarfile from the GMAO extranet server
        for colname in files:
            #-- sync Spire-GNSS files with GMAO server
            REMOTE = [*PATH,'grazing_angle_L2',colname]
            kwds = dict(DIRECTORY=DIRECTORY, TIMEOUT=TIMEOUT,
                MODE=MODE)
            output = multiprocess_sync(REMOTE, **kwds)
            #-- print the output string
            logging.info(output)
    else:
        #-- set multiprocessing start method
        ctx = mp.get_context("fork")
        #-- sync in parallel with multiprocessing Pool
        pool = ctx.Pool(processes=PROCESSES)
        #-- sync each Spire-GNSS data file
        out = []
        #-- get each tarfile from the GMAO extranet server
        for colname in files:
            #-- sync Spire-GNSS files with GMAO server
            REMOTE = [*PATH,'grazing_angle_L2',colname]
            kwds = dict(DIRECTORY=DIRECTORY, TIMEOUT=TIMEOUT,
                MODE=MODE)
            out.append(pool.apply_async(multiprocess_sync,
                args=(REMOTE,),kwds=kwds))
        #-- start multiprocessing jobs
        #-- close the pool
        #-- prevents more tasks from being submitted to the pool
        pool.close()
        #-- exit the completed processes
        pool.join()
        #-- print the output string
        for output in out:
            temp = output.get()
            logging.info(temp)

    #-- close log file and set permissions level to MODE
    if LOG:
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

#-- PURPOSE: wrapper for running the sync program in multiprocessing mode
def multiprocess_sync(*args, **kwds):
    try:
        output = http_pull_file(*args, **kwds)
    except Exception as e:
        #-- if there has been an error exception
        #-- print the type, value, and stack trace of the
        #-- current exception being handled
        logging.critical('process id {0:d} failed'.format(os.getpid()))
        logging.error(traceback.format_exc())
    else:
        return output

#-- PURPOSE: try extracting Spire files from a tar file
def http_pull_file(REMOTE, DIRECTORY=None, TIMEOUT=None, MODE=None):
    #-- regular expression pattern for extracting data from files
    rx = re.compile(r'(spire_gnss-r)_(L\d+)_(.*?)_(v\d+\.\d+)_'
        r'(\d{4})-(\d{2})-(\d{2})T(\d{2})-(\d{2})-(\d{2})_(.*?)\.nc$')
    #-- get BytesIO object containing data from tar file
    buffer = spire_toolkit.utilities.from_http(REMOTE,
        timeout=TIMEOUT, context=None)
    #-- open the monthly tar file
    tar1 = tarfile.open(fileobj=buffer, mode='r')
    mem1 = [f for f in tar1.getmembers() if f.name.endswith('tar')]
    #-- for each file within the tarfile
    output = ''
    for member in mem1:
        #-- print tarfile name to log
        output += '\t{0}\n'.format(member.name)
        #-- open the daily tarfile
        fileID = io.BytesIO(tar1.extractfile(member).read())
        tar2 = tarfile.open(fileobj=fileID, mode='r')
        mem2 = [f for f in tar2.getmembers() if f.name.endswith('nc')]
        #-- extract netCDF4 files to local directory
        for nc in mem2:
            #-- extract parameters from netCDF4 file
            MS,LV,PRD,VERS,YY,MM,DD,HH,MN,SS,AUX = rx.findall(nc.name).pop()
            SUBDIRECTORY = '{0}.{1}.{2}'.format(YY,MM,DD)
            #-- create output local directory if non-existent
            if not os.access(os.path.join(DIRECTORY,SUBDIRECTORY), os.F_OK):
                os.makedirs(os.path.join(DIRECTORY,SUBDIRECTORY), mode=MODE)
            #-- extract file
            tar2.extract(nc, path=os.path.join(DIRECTORY,SUBDIRECTORY))
            local_file = os.path.join(DIRECTORY,SUBDIRECTORY,nc.name)
            output += '\t\t{0} -->\n\t\t\t{1}\n'.format(nc.name,local_file)
            #-- use original modification time from tar file
            os.utime(local_file,(os.stat(local_file).st_atime,nc.mtime))
            #-- change the permissions mode
            os.chmod(local_file,MODE)
    #-- return the output string
    return output

#-- Main program that calls gmao_spire_gnss_sync()
def main():
   #-- Read the system arguments listed after the program
    parser = argparse.ArgumentParser(
        description="""Syncs Spire GNSS grazing angle altimetry data from
            the NASA Global Modeling and Assimilation Office (GMAO)
            """
    )
    #-- command line parameters
    parser.add_argument('--user','-U',
        type=str, default=os.environ.get('GMAO_USERNAME'),
        help='Username for NASA GMAO Extranet Login')
    parser.add_argument('--password','-W',
        type=str, default=os.environ.get('GMAO_PASSWORD'),
        help='Password for NASA GMAO Extranet Login')
    parser.add_argument('--netrc','-N',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.path.join(os.path.expanduser('~'),'.netrc'),
        help='Path to .netrc file for authentication')
    #-- working data directory
    parser.add_argument('--directory','-D',
        type=lambda p: os.path.abspath(os.path.expanduser(p)),
        default=os.getcwd(),
        help='Working data directory')
    #-- Spire data products to sync
    products = ['grzAlt','grzIce']
    parser.add_argument('--product','-p',
        type=int, nargs='+', choices=products, default=products,
        help='Spire data products to sync')
    #-- years of Spire data to sync
    now = time.gmtime()
    parser.add_argument('--year','-Y',
        type=int, nargs='+', default=range(2020,now.tm_year+1),
        help='Years of Spire data to sync')
    #-- run sync in series if processes is 0
    parser.add_argument('--np','-P',
        metavar='PROCESSES', type=int, default=0,
        help='Number of processes to use in file downloads')
    #-- connection timeout
    parser.add_argument('--timeout','-t',
        type=int, default=360,
        help='Timeout in seconds for blocking operations')
    #-- Output log file in form
    #-- GMAO_Spire_GNSS_sync_2002-04-01.log
    parser.add_argument('--log','-l',
        default=False, action='store_true',
        help='Output log file')
    #-- permissions mode of the directories and files synced (number in octal)
    parser.add_argument('--mode','-M',
        type=lambda x: int(x,base=8), default=0o775,
        help='Permission mode of directories and files synced')
    args,_ = parser.parse_known_args()

    #-- NASA GMAO Extranet hostname
    URS = 'gmao.gsfc.nasa.gov'
    #-- get NASA Earthdata credentials
    try:
        args.user,_,args.password = netrc.netrc(args.netrc).authenticators(URS)
    except:
        #-- check that NASA Earthdata credentials were entered
        if not args.user:
            prompt = 'Username for {0}: '.format(URS)
            args.user = builtins.input(prompt)
        #-- enter password securely from command-line
        if not args.password:
            prompt = 'Password for {0}@{1}: '.format(args.user,URS)
            args.password = getpass.getpass(prompt)

    #-- build an urllib opener for NASA GMAO Extranet
    #-- Add the username and password for NASA GMAO Extranet Login system
    opener = spire_toolkit.utilities.build_opener(args.user,args.password,
        password_manager=False,
        authorization_header=True,
        get_ca_certs=False,
        redirect=False)

    #-- post credentials to login to retrieve cookies
    LOGIN = posixpath.join('https://gmao.gsfc.nasa.gov','extranet','index.php')
    data = spire_toolkit.utilities.urlencode({'un':args.user, 'pw':args.password})
    request = spire_toolkit.utilities.urllib2.Request(LOGIN)
    response = spire_toolkit.utilities.urllib2.urlopen(request,
        data=data.encode('utf-8'))
    #-- verify url and cookies
    assert response.url
    assert opener.handlers[7].cookiejar

    #-- check internet connection before attempting to run program
    if spire_toolkit.utilities.check_connection(LOGIN):
        gmao_spire_gnss_sync(args.directory,
            PRODUCT=args.product,
            YEAR=args.year,
            PROCESSES=args.np,
            TIMEOUT=args.timeout,
            LOG=args.log,
            MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
