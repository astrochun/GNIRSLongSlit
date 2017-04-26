"""
reduce
======

Set of PyRAF routines to reduce GNIRS longslit spectra
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits
import glob # + on 26/04/2017

import numpy as np

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log # Mod on 26/04/2017

from pyraf import iraf

iraf.gemini(_doprint=0)
iraf.gemini.gnirs(_doprint=0)

log.info("Unlearning tasks")
iraf.gemini.unlearn()
iraf.gemini.gemtools.unlearn()
iraf.gemini.gnirs.unlearn()

iraf.gemini.nsheaders("gnirs")

# Set the display
iraf.set(stdimage="imt4096")

yes, no = 'yes', 'no' # + on 26/04/2017

def run(rawdir, bpm="gnirs$data/gnirsn_2012dec05_bpm.fits",
        silent=False, verbose=True):

    '''
    Main function to run the IRAF Gemini reduction package on GNIRS data

    Parameters
    ----------
    rawdir : str
      Path to raw files. Must end in a '/'

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 21 March 2017
    Modified by Chun Ly, 26 April 2017
     - Add code to run nsprepare
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    log.info("## Raw data is located in %s" % rawdir)

    # + on 26/04/2017
    all_lis = np.loadtxt(rawdir+'all.lis', dtype=type(str))
    n_all   = len(all_lis)
    
    # Step 1 | + on 26/04/2017
    log.info("## Preparing GNIRS data")

    nc_files = glob.glob(rawdir+'ncN*fits')
    n_nc     =  len(nc_files)

    if n_nc == n_all:
        log.warn("## File exists. Will not run nsprepare")
    else:
        fl_forcewcs = yes
        if n_nc == 0:
            iraf.gnirs.nsprepare(inimages="c@all.lis", rawpath=rawdir,
                                 bpm=bpm, shiftx=INDEF, shifty=INDEF,
                                 fl_forcewcs=fl_forcewcs)
        else:
            print 'blah'
    #endelse

    if silent == False: log.info('### End run : '+systime())
#enddef

