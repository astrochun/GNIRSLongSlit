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
import file_handling # + on 07/05/2017
import QA_wave_cal # + on 25/05/2017

co_filename = __file__ # + on 05/05/2017

yes, no = 'yes', 'no' # + on 26/04/2017

def run(rawdir, bpm="gnirs$data/gnirsn_2012dec05_bpm.fits",
        do_all=0, prepare=0, do_flat=0, do_arcs=0, wave_cal=0, skysub=0,
        fitcoords=0, combine=0, extract=0, silent=False, verbose=True):

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
     - Execution within rawdir rather than in parent directory
     - Use iraf_get_subset.check_prefix to check flats_rev files
    Modified by Chun Ly, 7 May 2017
     - Fix issue with call to nsprepare (seems to be a directory issue).
       Error given is:
        Warning: Cannot open file (tmpext43009_433)
        ERROR - GEMHEDIT: Image ncN20170218S0509[VAR] not found
        ERROR - GEMHEDIT: Image ncN20170218S0509[DQ] not found
       Problem seems to stem from the tmpext file being placed in the rawdir
       but iraf thinks it's in the directory before os.chdir()
       Reverting back so no os.chdir()
     - Modify call to nswavelength to get things to work within
       cdir (no need to change directory. Require copying files around
    Modified by Chun Ly, 15 May 2017
     - Use file_handling.mv_files instead of cp_files()
     - Simplify arc_list to include full path
     - Changes to call of file_handling.rm_files
    Modified by Chun Ly, 16 May 2017
     - Full path for flatfile, minor fixes
     - Add code to run skysub on telluric data with nsreduce (step5a)
     - Add code to run skysub on science data with nsreduce (step5b)
     - Add optional input keywords to run separate steps

     - Minor bug fix for tot0
     - Define all input list filenames at the beginning of this function
     - Switch input list filenames to include rawdir
     - Remove rawdir in multiple lines calling input list
     - Modify all calls to iraf_get_subset.check_prefix()
    Modified by Chun Ly, 17 May 2017
     - No need to run nsreduce on sky.lis since obj.lis include all frames
     - Added do_all keyword to simplify things. This will run all steps
     - Added additional log.info messages for each step
     - Add fitcoords keyword and code to execute nsfitcoords and nstransform

     - Add combine keyword and code to execute nscombine
     - Check if combined FITS file exists first before running nscombine

     - Add extract keyword and code to execute nsextract on telluric data

    Modified by Chun Ly, 20 May 2017
     - Use iraf.chdir to move around for nswavelength (step4)
    Modified by Chun Ly, 25 May 2017
     - Call QA_wave_cal.arc_check() to illustrate wavelength calibration
       on arc data
     - Call QA_wave_cal.OH_check() to illustrate wavelength calibration
       on object data
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    # + on 16/05/2017
    tot0 = sum([do_all, prepare, do_flat, do_arcs, wave_cal, skysub,
                fitcoords, combine, extract])
    if tot0 == 0:
        log.warn('## No GNIRS functions are being called!!!')
        log.warn('## Run with do_all=1 or either of these keywords set to 1 : ')
        log.warn('## prepare  : Execute nsprepare on all.lis')
        log.warn('## do_flat  : Create superflat')
        log.warn('## do_arcs  : nsreduce arc data')
        log.warn('## wave_cal : Wavelength calibration')
        log.warn('## skysub   : Skysubtraction - telluric and science data')
        log.warn('## fitcoords: nsfitcoords, nstransform on telluric and science data')
        log.warn('## combine  : nscombine telluric and science data')
        log.warn('## extract  : nsextract telluric data')
        return
    if do_all:
        prepare, do_flat, do_arcs = 1, 1, 1
        wave_cals, skysub = 1, 1
        extract, combine, fitcoords = 1, 1, 1

    cdir = os.getcwd()+'/' # + on 06/05/2017

    iraf.gemini.nsheaders("gnirs")

    # + on 26/04/2017
    timestamp = systime().replace(':','.')
    logfile   = rawdir+'gnirs_'+timestamp+'.log'
    iraf.gemini.gnirs.logfile = logfile

    log.info("## Raw data is located in : %s" % rawdir)

    log.info("## GNIRS logfile : "+logfile) # + on 05/05/2017

    # Save reduce.py for each run | + on 05/05/2017
    reduce_file = 'reduce_'+timestamp+'.py'
    log.info("## GNIRSLongSlit.reduce script : " + reduce_file)
    os.system('cp -a '+co_filename+' '+rawdir+reduce_file)

    # Check for cleanir files first | Later + on 26/04/2017
    c_files = glob.glob('cN*fits') # Mod on 06/05/2017
    if len(c_files) == 0:
        log.warn("## No cleanir files (cN*fits) available")
        log.warn("## Need to execute symlink.run()") # + on 05/05/2017
        log.warn("## ABORTING!!!")
        return
    #endif

    # Moved up on 16/05/2017
    all_list  = rawdir+'all.lis'
    flat_list = rawdir+'flat.lis'
    flats_rev = rawdir+'flat_rev.lis'
    arc_list  = rawdir+'arc.lis'
    tell_list = rawdir+'telluric.lis'
    obj_list  = rawdir+'obj.lis'
    sky_list  = rawdir+'sky.lis'
    flatfile  = rawdir+'final_flat.fits'

    # + on 17/05/2017
    tell_comb = rawdir+'tell_comb.fits'
    obj_comb  = rawdir+'obj_comb.fits'

    # + on 26/04/2017
    all_lis = np.loadtxt(all_list, dtype=type(str)) # Mod on 06/05/2017
    n_all   = len(all_lis)
    
    # Step 1 - Prepare GNIRS data | + on 26/04/2017
    if prepare:
        log.info("## Preparing GNIRS data")

        nc_files = glob.glob(rawdir+'ncN*fits') # Mod on 06/05/2017
        n_nc     =  len(nc_files)

        if n_nc == n_all:
            log.warn("## Files exist! Will not run nsprepare!!")
        else:
            fl_forcewcs = yes
            if n_nc == 0:
                # Mod on 26/04/2017 - Must specify full path
                # Mod on 06/05/2017
                inimages  = "c@"+all_list
                outimages = rawdir+"nc@"+all_list
                iraf.gnirs.nsprepare(inimages=inimages, rawpath=rawdir,
                                     outprefix='', outimages=outimages,
                                     bpm=bpm, shiftx="INDEF", shifty="INDEF",
                                     fl_forcewcs=fl_forcewcs, fl_correct=no)
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
    #end prepare

    # Step 2 - Create flat | + on 26/04/2017
    if do_flat:
        log.info("## Creating super flat")
        flats     = np.loadtxt(flat_list, dtype=type(str)) # Mod on 06/05/2017
        tmpflat   = rawdir+'tmpflat'
        if not exists(tmpflat):
            flat_files = ['nc'+file0+'[SCI,1]' for file0 in flats] # Mod on 06/05/2017
            if silent == False: log.info('## Writing : '+tmpflat)
            np.savetxt(tmpflat, flat_files, fmt='%s')
        else:
            log.warn('## File exists!!! : '+tmpflat)

        # Compute stats on flat; Use reliable ones + on 04/05/2017
        flat_sig = iraf.imstatistic(images=rawdir+'@'+tmpflat,
                                    lower=0, nclip=5,format='no',
                                    fields="image,npix,mean,stddev,min,max",
                                    Stderr=1)

        npix = [int(str0.split('  ')[1])   for str0 in flat_sig]
        mean = [float(str0.split('  ')[2]) for str0 in flat_sig]
        std  = [float(str0.split('  ')[3]) for str0 in flat_sig]

        avg  = np.average(mean)
        rat0 = 100.0*np.absolute(1-mean/avg)
        good = np.where(rat0 <= 0.5) # If difference is more than 0.5%
        if len(good) > 0:
            log.info('## Flat files to use : ')
            log.info(', '.join(flats[good]))
            np.savetxt(flats_rev, flats[good], fmt='%s')

        # + on 04/05/2017 | Mod on 05/05/2017
        # Mod on 06/05/2017
        do_run = iraf_get_subset.check_prefix('rnc', flats_rev)
        if do_run:
            # Mod on 06/05/2017
            iraf.gnirs.nsreduce(rawdir+'nc@'+flats_rev, outprefix='',
                                outimages=rawdir+'rnc@'+flats_rev,
                                fl_sky=no, fl_cut=yes, fl_flat=no, fl_dark=no,
                                fl_nsappwave=no)
        else:
            log.warn("## Files exist! Will not run nsreduce!!")

        # + on 05/05/2017
        if not exists(flatfile):
            iraf.gnirs.nsflat(rawdir+'rnc@'+flats_rev, flatfile=flatfile) # Mod on 06/05/2017
        else:
            log.warn('## File exists!!! : '+flatfile)
            log.warn('## Will not run nsflat')
    #end do_flat

    # Step 3 : Reduce arcs | + on 05/05/2017
    arcs = np.loadtxt(arc_list, dtype=type(str)) # Mod on 06/05/2017

    if do_arcs:
        log.info("## Reducing arc data")
        do_run = iraf_get_subset.check_prefix('rnc', arc_list)
        if do_run:
            # Mod on 06/05/2017
            iraf.gnirs.nsreduce(rawdir+'nc@'+arc_list, outprefix='',
                                outimages=rawdir+'rnc@'+arc_list,
                                fl_sky=no, fl_cut=yes, fl_flat=no,
                                fl_dark=no) #fl_nsappwave=no)
        else:
            log.warn('## File exists!!!')
            log.warn('## Will not run nsreduce on arc data')
    #end do_arcs

    # Step 4 : Perform wavelength calibration | + on 05/05/2017
    if wave_cal:
        iraf.chdir(rawdir) # + on 20/05/2017
        log.info("## Performing non-interactive wavelength calibration on arc data")
        do_run = iraf_get_subset.check_prefix('wrnc', arc_list)
        if do_run:
            # Mod on 06/05/2017
            iraf.gnirs.nswavelength('rnc@'+arc_list, outprefix='',
                                    outspectra='wrnc@'+arc_list,
                                    coordlist="gnirs$data/lowresargon.dat",
                                    database='database/',
                                    fl_inter=no, cradius=20, threshold=50.0,
                                    order=2)
        else:
            log.warn('## Files exist!!!')
            log.warn('## Will not run nswavelength on rnc arc data')
        iraf.chdir(cdir)

        QA_wave_cal.arc_check(rawdir, arcs=arcs) # + on 25/05/2017

    # Step 5a : Sky subtract telluric data | + on 16/05/2017
    if skysub:
        log.info("## Performing sky subtraction on telluric data")
        do_run = iraf_get_subset.check_prefix('rnc', tell_list)
        if do_run:
            iraf.gnirs.nsreduce(rawdir+'nc@'+tell_list, outprefix='',
                                outimages=rawdir+'rnc@'+tell_list,
                                fl_nsappwave=no, fl_sky=yes, fl_flat=yes,
                                flatimage=flatfile)
        else:
            log.warn('## Files exist!!!')
            log.warn('## Will not run nsreduce on nc telluric data')

        # Step 5b : Sky subtract science data | + on 16/05/2017
        log.info("## Performing sky subtraction on science data")
        do_run = iraf_get_subset.check_prefix('rnc', obj_list)
        if do_run:
            iraf.gnirs.nsreduce(rawdir+'nc@'+obj_list, outprefix='',
                                outimages=rawdir+'rnc@'+obj_list,
                                fl_nsappwave=no, fl_sky=yes,
                                skyimages=rawdir+'nc@'+sky_list,
                                fl_flat=yes, flatimage=flatfile)
        else:
            log.warn('## Files exist!!!')
            log.warn('## Will not run nsreduce on nc sci data')
    #end skysub

    # Step 6a : Apply wavelength solution to telluric data | + on 17/05/2017
    if fitcoords:
        # Telluric data
        iraf.chdir(rawdir)
        do_run = iraf_get_subset.check_prefix('frnc', tell_list)
        if do_run:
            log.info("## Running nsfitcoords on telluric data")
            iraf.gnirs.nsfitcoords('rnc@'+tell_list, outprefix='',
                                   outspectra='frnc@'+tell_list,
                                   lamp='wrnc'+arcs[0],
                                   database='database/')
        else:
            log.warn('## Files exist!!!')
            log.warn('## Will not run nsfitcoords on rnc telluric data')

        do_run = iraf_get_subset.check_prefix('tfrnc', tell_list)
        if do_run:
            log.info("## Running nstransform on telluric data")
            iraf.gnirs.nstransform('frnc@'+tell_list, outprefix='',
                                   outspectra='tfrnc@'+tell_list,
                                   database='database/')
        else:
            log.warn('## Files exist!!!')
            log.warn('## Will not run nstransform on frnc telluric data')

        # Step 6b : Apply wavelength solution to science data | + on 17/05/2017
        # Science data
        do_run = iraf_get_subset.check_prefix('frnc', obj_list)
        if do_run:
            log.info("## Running nsfitcoords on science data")
            iraf.gnirs.nsfitcoords('rnc@'+obj_list, outprefix='',
                                   outspectra='frnc@'+obj_list,
                                   lamp='wrnc'+arcs[0],
                                   database='database/')
        else:
            log.warn('## Files exist!!!')
            log.warn('## Will not run nsfitcoords on rnc science data')

        do_run = iraf_get_subset.check_prefix('tfrnc', obj_list)
        if do_run:
            log.info("## Running nstransform on science data")
            iraf.gnirs.nstransform('frnc@'+obj_list, outprefix='',
                                   outspectra='tfrnc@'+obj_list,
                                   database='database/')
        else:
            log.warn('## Files exist!!!')
            log.warn('## Will not run nstransform on frnc science data')

        iraf.chdir(cdir)
        QA_wave_cal.OH_check(rawdir) # + on 25/05/2017
    #end fitcoords

    # Step 7: Combine 2-D spectra | + on 17/05/2017
    if combine:
        if not exists(tell_comb):
            log.info("## Running nscombine on telluric data")
            iraf.gnirs.nscombine(rawdir+'tfrnc@'+tell_list, output=tell_comb,
                                 fl_cross=yes, tolerance=0.1)
        else:
            log.warn('## File exists : '+tell_comb+' !!!')
            log.warn('## Will not run nscombine on tfrnc telluric data')

        if not exists(obj_comb):
            log.info("## Running nscombine on science data")
            iraf.gnirs.nscombine(rawdir+'tfrnc@'+obj_list, output=obj_comb,
                                 fl_cross=yes, tolerance=0.1)
        else:
            log.warn('## File exists : '+obj_comb+' !!!')
            log.warn('## Will not run nscombine on tfrnc science data')

    # Step 8: Extract 1-D spectra | + on 17/05/2017
    if extract:
        outspec = tell_comb.replace('tell','xtell')
        if not exists(outspec):
            log.info("## Running nsextract on telluric data")
            iraf.gnirs.nsextract(tell_comb, outspectra=outspec,
                                 database=rawdir+'database/')
        else:
            log.warn('## File exists : '+outspec+' !!!')
            log.warn('## Will not run nsextract on comb_tell.fits')

    #os.chdir(cdir) # + on 06/05/2017

    if silent == False: log.info('### End run : '+systime())
#enddef

