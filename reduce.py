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
     - Add code to run nsprepare (step1)
     - Check for cN*fits files (cleanir or symlink files)
     - Load in iraf_get_subset package
     - Fix path bug when calling nsprepare
     - Add warning if not all nsprepare files are available
     - Define logfile for GNIRS
     - Add code to compute statistics (step2)
    Modified by Chun Ly, 4 May 2017
     - Compute statistics on flats and use the most reliable ones
       (exclude outliers)
     - Run nsreduce on flats
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    iraf.gemini.nsheaders("gnirs")

    # + on 26/04/2017
    timestamp = systime()
    iraf.gemini.gnirs.logfile = rawdir+'gnirs_'+timestamp+'.log'

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
            # Mod on 26/04/2017 - Must specify full path
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

    # Step 2 - Create flat | + on 26/04/2017
    flats = np.loadtxt(rawdir+'flat.lis', dtype=type(str))
    tmpflat = rawdir+'tmpflat'
    if not exists(tmpflat):
        flat_files = [rawdir+'nc'+file0+'[SCI,1]' for file0 in flats]
        if silent == False: log.info('## Writing : '+tmpflat)
        np.savetxt(tmpflat, flat_files, fmt='%s')
    else:
        if silent == False: log.warn('## File exists!!! : '+tmpflat)

    # Compute stats on flat; Use reliable ones + on 04/05/2017
    flat_sig = iraf.imstatistic(images='@'+tmpflat, lower=0, nclip=5,
                                fields="image,npix,mean,stddev,min,max",
                                format='no', Stderr=1)
    npix = [int(str0.split('  ')[1])   for str0 in flat_sig]
    mean = [float(str0.split('  ')[2]) for str0 in flat_sig]
    std  = [float(str0.split('  ')[3]) for str0 in flat_sig]

    avg  = np.average(mean)
    rat0 = 100.0*np.absolute(1-mean/avg)
    good = np.where(rat0 <= 0.5)
    flats_rev = rawdir+'flat_rev.lis'
    if len(good) > 0:
        log.info('## Flat files to use : ')
        log.info('\n'.join(flats[good]))
        np.savetxt(flats_rev, flats[good], fmt='%s')

    # + on 04/05/2017
    iraf.gnirs.nsreduce(rawdir+'nc@'+flats_rev,
                        outimages=rawdir+'rnc@'+flats_rev, outprefix='',
                        fl_sky=no, fl_cut=yes, fl_flat=no, fl_dark=no,
                        fl_nsappwave=no)

    if silent == False: log.info('### End run : '+systime())
#enddef

