"""
symlink
=======

Create symbolic links between raw data for prefix handling of cleanir products

Code checks to see if a cNYYYYMMDDSNNNN.fits exists. If not, it will provide
do a symbolic link
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy.table import Table
from astropy import log

from . import gnirs_2017a

import dir_check # + on 23/03/2017

import glog # + on 08/01/2018

def get_files(path0, silent=False, verbose=True):
    '''
    Simple function to get names of raw files

    Parameters
    ----------
    path0 : str
     Path to raw FITS file. Must include '/' at the end

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    files : list
      List containing the filename for FITS raw files

    n_files : int
      Size of [files]

    Notes
    -----
    Created by Chun Ly, 22 March 2017
    Modified by Chun Ly, 5 June 2017
     - Fix minor bug: infile -> infile0
    Modified by Chun Ly, 8 January 2018
     - Import glog and call for stdout and ASCII logging
    '''

    logfile  = path0+'symlink.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin get_files : '+systime())

    infile0 = path0+'all.lis'
    if not exists(infile0):
        mylogger.warn('File does not exists!!! : '+infile0)
        mylogger.warn('EXITING!!!')
        return

    if silent == False:
        mylogger.info('Reading : '+infile0)
    files   = np.loadtxt(infile0, dtype=type(str)).tolist()
    n_files = len(files)

    if silent == False: mylogger.info('### End get_files : '+systime())

    return files, n_files
#enddef



def delete(path0, silent=False, verbose=True):
    '''
    Remove all symbolic links in given path containing the raw GNIRS data

    Parameters
    ----------
    path0 : str
     Path to raw FITS file. Must include '/' at the end

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 22 March 2017
    Modified by Chun Ly, 8 January 2018
     - Import glog and call for stdout and ASCII logging
    '''

    # + on 08/01/2018
    logfile  = path0+'symlink.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin delete : '+systime())

    files, n_files = get_files(path0, silent=silent, verbose=verbose)

    for nn in xrange(n_files):
        c_file = path0+'c'+files[nn]
        if exists(c_file):
            if os.path.islink(c_file) == True:
                cmd0 = 'rm '+c_file
                log.info(cmd0)
                os.system(cmd0)
            else:
                mylogger.info('File is from cleanir: '+c_file)

    if silent == False: mylogger.info('### End delete : '+systime())
#enddef

def run(path0, silent=False, verbose=True):
    '''
    Main function to run for a given path containing the raw GNIRS data

    Parameters
    ----------
    path0 : str
     Path to raw FITS file. Must include '/' at the end

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 22 March 2017
    Modified by Chun Ly, 23 March 2017
     - Call dir_check.main() to handle multiple date directories
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    # + on 23/03/2017
    dir_list, list_path = dir_check.main(path0, silent=silent, verbose=verbose)

    # Mod on 23/03/2017
    for path in list_path:
        files, n_files = get_files(path, silent=silent, verbose=verbose)

        for nn in xrange(n_files):
            c_file = path+'c'+files[nn]
            if exists(c_file):
                log.warn('File exists : '+c_file)
            else:
                cmd0 = 'ln -fs '+files[nn]+' '+c_file
                if silent == False:
                    log.info(cmd0)
                os.system(cmd0)
        #endfor
    #endfor

    if silent == False: log.info('### End run : '+systime())
#enddef

def zcalbase_gal_gemini_2017a():
    '''
    Function to execute run() on each set of GNIRS 2017A observation set
    to generate symbolic links for raw data that do not need cleanir executed

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 22 March 2017
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = gnirs_2017a

    for target in targets0:
        run(path0+target+'/')

#enddef
