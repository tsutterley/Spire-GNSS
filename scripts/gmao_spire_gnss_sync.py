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
    -P X, --password X: password for NASA GMAO Extranet Login
    -N X, --netrc X: path to .netrc file for authentication
    -D X, --directory X: working data directory
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
import tarfile
import builtins
import argparse
import posixpath
import spire_toolkit.utilities

#-- PURPOSE: Syncs Spire GNSS grazing angle altimetry data
def gmao_spire_gnss_sync(DIRECTORY, TIMEOUT=None, LOG=False, MODE=0o775):
    #-- create log file with list of synchronized files (or print to terminal)
    if LOG:
        #-- format: GMAO_Spire_GNSS_sync_2002-04-01.log
        today = time.strftime('%Y-%m-%d',time.localtime())
        LOGFILE = 'GMAO_Spire_GNSS_sync_{0}.log'.format(today)
        fid1 = open(os.path.join(DIRECTORY,LOGFILE),'w')
        print('GMAO Spire GNSS Sync Log ({0})'.format(today), file=fid1)
    else:
        #-- standard output (terminal output)
        fid1 = sys.stdout

    #-- remote host for Spire GNSS data
    HOST = 'https://gmao.gsfc.nasa.gov'
    #-- open connection with GMAO extranet server at remote directory
    PATH = [HOST,'extranet','collab','spire_team']
    files = spire_toolkit.utilities.gmao_list(PATH, timeout=TIMEOUT,
        adddirlink="grazing_angle_L2", sort=True)
    #-- regular expression pattern for extracting data from files
    rx = re.compile(r'(spire_gnss-r)_(L\d+)_(.*?)_(v\d+\.\d+)_'
        r'(\d{4})-(\d{2})-(\d{2})T(\d{2})-(\d{2})-(\d{2})_(.*?)\.nc$')
    #-- get each tarfile from the GMAO extranet server
    for colname in files:
        #-- print tarfile name to log
        print(colname, file=fid1)
        PATH = [*PATH,'grazing_angle_L2',colname]
        buffer = spire_toolkit.utilities.from_http(PATH,
            timeout=TIMEOUT, context=None)
        #-- open the monthly tar file
        tar1 = tarfile.open(fileobj=buffer, mode='r')
        mem1 = [f for f in tar1.getmembers() if f.name.endswith('tar')]
        #-- for each file within the tarfile
        for member in mem1:
            #-- print tarfile name to log
            print('\t{0}'.format(member.name), file=fid1)
            #-- open the daily tarfile
            fileID = io.BytesIO(tar1.extractfile(member).read())
            tar2 = tarfile.open(fileobj=fileID, mode='r')
            mem2 = [f for f in tar2.getmembers() if f.name.endswith('nc')]
            #-- extract netCDF4 files to local directory
            for nc in mem2:
                print('\t\t{0} -->'.format(nc.name), file=fid1)
                #-- extract parameters from netCDF4 file
                MS,LV,PRD,VERS,YY,MM,DD,HH,MN,SS,AUX = rx.findall(nc.name).pop()
                SUBDIRECTORY = '{0}.{1}.{2}'.format(YY,MM,DD)
                if not os.access(os.path.join(DIRECTORY,SUBDIRECTORY), os.F_OK):
                    os.makedirs(os.path.join(DIRECTORY,SUBDIRECTORY), mode=MODE)
                #-- extract file
                tar2.extract(nc, path=os.path.join(DIRECTORY,SUBDIRECTORY))
                local_file = os.path.join(DIRECTORY,SUBDIRECTORY,nc.name)
                print('\t\t\t{0} -->'.format(local_file), file=fid1)
                #-- change the permissions mode
                os.chmod(local_file,MODE)

    #-- close log file and set permissions level to MODE
    if LOG:
        fid1.close()
        os.chmod(os.path.join(DIRECTORY,LOGFILE), MODE)

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
        gmao_spire_gnss_sync(args.directory, TIMEOUT=args.timeout,
            LOG=args.log, MODE=args.mode)

#-- run main program
if __name__ == '__main__':
    main()
