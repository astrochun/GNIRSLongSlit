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

# Set the display
iraf.set(stdimage="imt4096")

import iraf_get_subset # + on 26/04/2017

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
     - Check for cN*fits files (cleanir or symlink files)
     - Load in iraf_get_subset package
     - Fix path bug when calling nsprepare
     - Add warning if not all nsprepare files are available
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    iraf.gemini.nsheaders("gnirs")

    log.info("## Raw data is located in %s" % rawdir)

    # Check for cleanir files first | Later + on 26/04/2017
    c_files = glob.glob(rawdir+'cN*fits')
    if len(c_files) == 0:
        log.warn("## No cleanir files (cN*fits) available")
        log.warn("## ABORTING!!!")
        return
    #endif

    # + on 26/04/2017
    all_lis = np.loadtxt(rawdir+'all.lis', dtype=type(str))
    n_all   = len(all_lis)
    
    # Step 1 - Prepare GNIRS data | + on 26/04/2017
    log.info("## Preparing GNIRS data")

    nc_files = glob.glob(rawdir+'ncN*fits')
    n_nc     =  len(nc_files)

    if n_nc == n_all:
        log.warn("## Files exist! Will not run nsprepare!!")
    else:
        fl_forcewcs = yes
        if n_nc == 0:
            # Mod on 26/04/2017 - Most specify full path
            inimages  = "c@"+rawdir+"all.lis"
            outimages = "nc@"+rawdir+"all.lis"
            iraf.gnirs.nsprepare(inimages=inimages, rawpath=rawdir,
                                 outimages=outimages, bpm=bpm,
                                 shiftx="INDEF", shifty="INDEF",
                                 fl_forcewcs=fl_forcewcs)
        else:
            '''
            Warns if files do not exist. Need to implement a way to run
            nsprepare for a subset of data (need to include certain frames)
            '''
            log.warn("## The following files do not exist: ")
            iraf_get_subset.main(rawdir, 'nc', all_lis=all_lis, silent=True)
            #outfile = 'nc_sub.lis'
            #iraf_get_subset.main(rawdir, 'nc', outfile, all_lis=all_lis)
            #iraf.gnirs.nsprepare(inimages="c@"+outfile, rawpath=rawdir,
            #                     bpm=bpm, shiftx="INDEF", shifty="INDEF",
            #                     fl_forcewcs=fl_forcewcs)
    #endelse

    if silent == False: log.info('### End run : '+systime())
#enddef

