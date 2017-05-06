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

co_filename = __file__ # + on 05/05/2017

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
    Modified by Chun Ly, 5 May 2017
     - Check for rnc (nsreduce) files for flats
     - Run nsflat on flats
     - Add code to run nsreduce on arcs (step3)
     - Save reduce.py for each execution of reduce.run()
    Modified by Chun Ly, 6 May 2017
     - Add code to run nswavelength on arcs (step4)
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    iraf.gemini.nsheaders("gnirs")

    # + on 26/04/2017
    timestamp = systime()
    logfile   = rawdir+'gnirs_'+timestamp+'.log'
    iraf.gemini.gnirs.logfile = logfile

    log.info("## Raw data is located in : %s" % rawdir)

    log.info("## GNIRS logfile : "+logfile) # + on 05/05/2017

    # Save reduce.py for each run | + on 05/05/2017
    reduce_file = rawdir+'reduce_'+timestamp+'.py'
    log.info("## GNIRSLongSlit.reduce script : " + reduce_file)
    os.system('cp -a '+co_filename+' '+reduce_file)

    # Check for cleanir files first | Later + on 26/04/2017
    c_files = glob.glob(rawdir+'cN*fits')
    if len(c_files) == 0:
        log.warn("## No cleanir files (cN*fits) available")
        log.warn("## Need to execute symlink.run()") # + on 05/05/2017
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
        log.warn('## File exists!!! : '+tmpflat)

    # Compute stats on flat; Use reliable ones + on 04/05/2017
    flat_sig = iraf.imstatistic(images='@'+tmpflat, lower=0, nclip=5,
                                fields="image,npix,mean,stddev,min,max",
                                format='no', Stderr=1)
    npix = [int(str0.split('  ')[1])   for str0 in flat_sig]
    mean = [float(str0.split('  ')[2]) for str0 in flat_sig]
    std  = [float(str0.split('  ')[3]) for str0 in flat_sig]

    avg  = np.average(mean)
    rat0 = 100.0*np.absolute(1-mean/avg)
    good = np.where(rat0 <= 0.5) # If difference is more than 0.5%
    flats_rev = rawdir+'flat_rev.lis'
    if len(good) > 0:
        log.info('## Flat files to use : ')
        log.info(', '.join(flats[good]))
        np.savetxt(flats_rev, flats[good], fmt='%s')

    # + on 04/05/2017 | Mod on 05/05/2017
    rnc_files = glob.glob(rawdir+'rncN*fits')
    n_rnc      =  len(rnc_files)
    if n_nc >= len(good):
        log.warn("## Files exist! Will not run nsreduce!!")
    else:
        iraf.gnirs.nsreduce(rawdir+'nc@'+flats_rev, outprefix='',
                            outimages=rawdir+'rnc@'+flats_rev, fl_sky=no,
                            fl_cut=yes, fl_flat=no, fl_dark=no, fl_nsappwave=no)

    # + on 05/05/2017
    flatfile = rawdir+'final_flat.fits'
    if not exists(flatfile):
        iraf.gnirs.nsflat(rawdir+'rnc@'+flats_rev, flatfile=flatfile)
    else:
        log.warn('## File exists!!! : '+flatfile)
        log.warn('## Will not run nsflat')


    # Step 3 : Reduce arcs | + on 05/05/2017
    arc_list = 'arc.lis'
    arcs     = np.loadtxt(rawdir+arc_list, dtype=type(str))

    do_run = iraf_get_subset.check_prefix(rawdir, 'rnc', arc_list)
    if do_run:
        iraf.gnirs.nsreduce(rawdir+'nc@'+rawdir+arc_list, outprefix='',
                            outimages=rawdir+'rnc@'+rawdir+arc_list,
                            fl_sky=no, fl_cut=yes, fl_flat=no,
                            fl_dark=no) #fl_nsappwave=no)
    else:
        log.warn('## File exists!!! : '+flatfile)
        log.warn('## Will not run nsreduce on arc data')

    # Step 4 : Perform wavelength calibration | + on 05/05/2017
    do_run = iraf_get_subset.check_prefix(rawdir, 'wrnc', arc_list)
    if do_run:
        iraf.gnirs.nswavelength(rawdir+'rnc@'+rawdir+arc_list, outprefix='',
                                outspectra=rawdir+'wrnc@'+rawdir+arc_list,
                                coordlist="gnirs$data/lowresargon.dat",
                                database=rawdir+'database/',
                                fl_inter=no, cradius=20, threshold=50.0,
                                order=2)
    else:
        log.warn('## Files exist!!!')
        log.warn('## Will not run nswavelength on rnc arc data')

    if silent == False: log.info('### End run : '+systime())
#enddef

