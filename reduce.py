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
from pylab import subplots_adjust

from astropy.table import Table
from astropy import log # Mod on 26/04/2017

from astropy.stats import sigma_clipped_stats # + on 16/06/2017
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
# import file_handling   # + on 07/05/2017
import QA_wave_cal     # + on 25/05/2017
import OH_stack        # + on 13/07/2017
import wave_cal_script # + on 12/11/2017

# + on 14/09/2017
import remove_bias_level
import examine_median

# + on 20/09/2017
from check_path import main as check_path

import glog # + on 10/12/2017

# + on 05/05/2017
# Mod on 16/01/2018: Bug fix - copy ASCII py code instead of binary .pyc
co_filename = __file__.replace('.pyc','.py')
yes, no = 'yes', 'no' # + on 26/04/2017

n_sp_pix = 1022 # + on 12/06/2017

def compute_weights(rawdir, out_pdf='', silent=True, verbose=False):
    '''
    Uses alignment star to determine weights for each science frame

    Parameters
    ----------
    rawdir : str
      Path to raw files. Must end in a '/'

    out_pdf : str
      Filename for output PDF. Do NOT include full path
      Default: 'compute_weights.pdf'

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 3 June 2017
    Modified by Chun Ly, 14 September 2017
     - Change prefix to use the bias-subtracted frames - from remove_bias_level()
    Modified by Chun Ly, 20 September 2017
     - Call check_path()
    '''

    if silent == False: log.info('### Begin compute_weights : '+systime())

    rawdir = check_path(rawdir) # + on 20/09/2017

    out_pdf = rawdir+'compute_weights.pdf' if out_pdf == '' else \
              rawdir+out_pdf

    obj_list  = rawdir+'obj.lis'

    do_run = iraf_get_subset.check_prefix('e', obj_list, path=rawdir)
    if do_run:
        iraf.gnirs.nsextract(rawdir+'tfrbnc@'+obj_list,
                             outspectra=rawdir+'e@'+obj_list, nsum=50,
                             database=rawdir+'database/')

    e_files = glob.glob(rawdir+'eN*fits')
    n_files = len(e_files)

    weights0 = np.zeros(n_files)
    median0 = np.zeros(n_files)

    fig, ax = plt.subplots()

    for nn in xrange(n_files):
        print '## Reading : ', e_files[nn]
        t_hdu = fits.open(e_files[nn])
        spec1d = t_hdu['SCI'].data
        hdr    = t_hdu['SCI'].header

        crval1 = hdr['CRVAL1']
        cd1_1  = hdr['CD1_1']
        etime  = t_hdu[0].header['EXPTIME']

        wave = (crval1 + cd1_1*np.arange(hdr['NAXIS1']))/1E4 # Microns
        cts  = spec1d/etime

        median0[nn] = np.median(cts)
        ax.plot(wave, cts, label=os.path.basename(e_files[nn]))
    #endfor

    ax.set_xlabel(r'Wavelength ($\mu$m)', fontsize=16)
    ax.set_ylabel('Flux (cts/s)', fontsize=16)
    ax.set_ylim([0.0, max(median0)*1.25])
    ax.minorticks_on()

    subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99)
    ax.legend(loc='upper right', fontsize=12, frameon=False)
    fig.set_size_inches(11,8)
    fig.savefig(out_pdf)
    if silent == False: log.info('### End compute_weights : '+systime())
#enddef

def normalize_flat(flatfile, out_pdf='', rawdir='', mylogger=None, silent=True,
                   verbose=False):
    '''
    Plot up flat average and fix weird flats (broad dip)

    Parameters
    ----------
    flatfile : str
      Filename of final flat. Must provide full path

    silent : boolean
      Turns off stdout messages. Default: True

    verbose : boolean
      Turns on additional stdout messages. Default: False

    Returns
    -------
    PDF plot. Default: 'normalized_flat.pdf'

    Notes
    -----
    Created by Chun Ly, 16 June 2017

    Modified by Chun Ly, 10 December 2017
     - Implement glog logging, allow rawdir and mylogger keywords
    '''

    # + on 10/12/2017
    if type(mylogger) == type(None):
        logfile  = rawdir+'reduce.log'
        mylogger = glog.log0(logfile)._get_logger()

    if silent == False:
        mylogger.info('Begin normalize_flat : '+systime())

    path0 = os.path.dirname(flatfile)+'/'

    if silent == False: mylogger.info('Reading : '+flatfile)
    hdu0 = fits.open(flatfile)
    im0  = hdu0['SCI'].data

    fig, ax = plt.subplots()

    mean, median, std = sigma_clipped_stats(im0, sigma=2, iters=10, axis=0)
    ax.plot(mean,       'b--', label='mean')
    ax.plot(median,     'k--', label='median')
    ax.plot(median+std, 'g:', label=r'$\sigma$')
    ax.plot(median-std, 'g:')
    ax.set_xlabel('X [pixels]')
    ax.set_ylabel('Normalized value')
    ax.minorticks_on()
    ax.set_xlim([0,700])
    ax.legend(loc='upper right', frameon=False)
    subplots_adjust(left=0.1, right=0.975, top=0.99, bottom=0.1)
    if out_pdf == '': out_pdf = path0 + 'normalize_flat.pdf'

    if silent == False: mylogger.info('Writing : '+out_pdf)
    fig.savefig(out_pdf)

    if silent == False: mylogger.info('End normalize_flat : '+systime())

#enddef

def computeStatistics(rawdir, flat_files):
    '''
    rawdir : str
      Path to raw files. Must end in a '/'
    
    Parameters
    ---------
    flat_files : list(str)
        list of flat files to be analyzed

    Notes
    -----
    Modified by Chun Ly, 20 November 2017
     - Fix bug with not having full path for flat_files list
    '''

    stats = []
    for image in flat_files:
        hdu = fits.open(rawdir + image) # Mod on 20/11/2017
        data = hdu['SCI'].data[:,160:800] #160-800 pixel range to exclude bias level
        stats.append(sigma_clipped_stats(data))
    stats = np.array(stats)
    mean = stats[:,0]
    median = stats[:,1]
    stddev = stats[:,2]
    avg = np.average(mean)
    rat0 = 100.0*np.absolute(1-mean/avg)
    goodFlats = np.where(rat0 <= 0.5)
    return goodFlats

def run(rawdir, bpm="gnirs$data/gnirsn_2012dec05_bpm.fits",
        do_all=0, prepare=0, do_flat=0, do_arcs=0, wave_cal=0, skysub=0,
        fitcoords=0, combine=0, extract=0, tell_corr=0, calib_line='OH',
        silent=False, verbose=True):

    '''
    Main function to run the IRAF Gemini reduction package on GNIRS data

    Parameters
    ----------
    rawdir : str
      Path to raw files. Must end in a '/'

    calib_line : str
      Indication for using OH sky lines or arc lines calibration
      Either 'OH' or 'arc'. Default: 'OH'

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
    Modified by Chun Ly, 2 June 2017
     - Run nsreduce without skysub for wavelength check using OH skylines
     - Fix call to iraf_get_subset.check_prefix()
     - Fix bug: variable name incorrect
     - Run nsfitcoords and nstransform on un-skysubtracted science data
       for wavelength check using OH skylines
     - Call QA_wave_cal.OH_check to produce OH_check.raw.pdf with
       un-skysubtracted science data
    Modified by Chun Ly, 6 June 2017
     - Handle flat-fielding for X-band data (no flat applied)
    Modified by Chun Ly, 7 June 2017
     - No flatfielding for telluric X-band data (for consistency)
     - Fix minor bug: Check for cN*fits file in rawdir
    Modified by Chun Ly, 12 June 2017
     - Pass initial guesses for CRVAL,CDELT,CRPIX to nswavelength
    Modified by Chun Ly, 15 June 2017
     - Run nsfitcoords and nstransform on arc data to examine wavelength
       calibration
    Modified by Chun Ly, 16 June 2017
     - Call QA_wave_cal.arc_check2 to produce arc_check2.pdf with
       transformed arc data
    Modified by Chun Ly, 26 June 2017
     - Change sign for CDELT
    Modified by Chun Ly, 11 July 2017
     - Call normalize_flat()
     - Logging information
     - Minor bug: wave_cals -> wave_cal
    Modified by Chun Ly, 13 July 2017
     - Add call to OH_stack.run() in wave_cal to generate stacked OH data
    Modified by Chun Ly, 10 September 2017
     - Found source of problem for flatfielding.  Seems like nsreduce is
       multiplying science frames by the flat rather than dividing. Creating
       a flat that is the inverse from nsflat to use for flatfielding
     - Generate normalize_flat plot using nsflat's flat
     - Use inverted flat for nsreduce call of telluric and science data
       with fl_flat='yes' to fix flatfielding bug
    Modified by Chun Ly, 14 September 2017
     - Call remove_bias_level()
     - Flag CRs in nsprepare
     - Fix typo in call to remove_bias_level()
     - Change prefix to use the bias-subtracted frames - from remove_bias_level()
     - Remove bias level in flats; Use bias-subtracted flats for superflat
     - Remove previous code involving inverted flat
     - Call examine_median for orig and bias-subtracted cases
    Modified by Chun Ly, 15 September 2017
     - Apply flatfield to bnc OH images
     - Call examine_median for skysub and flat cases
    Modified by Chun Ly, 20 September 2017
     - Call check_path()
    Modified by Scott McKinley, 5 October 2017
     - Created computeStatistics()
     - removed call to iraf imstatistics with do_flat
    Modified by Chun Ly, 7 November 2017
     - Create and write average arcs stack
    Modified by Chun Ly, 8 November 2017
     - Bug with average arcs stack - Not formatted for gnirs pipeline
    Modified by Chun Ly, 11 November 2017
     - Modify wave_cal sub-routine to handle interactive fitting from PyRAF
     - Change call to QA_wave_cal.arc_check() to use stacked arc products
    Modified by Chun Ly, 12 November 2017
     - Use stacked arc products in call to nsfitcoords and nswavelength
     - Change call to QA_wave_cal.arc_check2() to include stacked arc products
     - Call wave_cal_script to generate .py script for interactive wavelength
       calibration with arc stacked data
    Modified by Chun Ly, 15 November 2017
     - Bug fix: 'arcs' -> 'arc'
     - Bug fix: 'arcs' -> 'arc' for filename

    Modified by Chun Ly, 16 November 2017
     - Minor documentation
    Modified by Chun Ly, 17 November 2017
     - Call wave_cal() and transform() of OH_stack
    Modified by Chun Ly, 20 November 2017
     - Pass rawdir into computeStatistics()
     - Always define flat_files list in do_flat if statement
     - Pass cdir into OH_stack.wave_cal()
    Modified by Chun Ly, 22 November 2017
     - Call OH_stack.plot_spec in wave_cal if statement
    Modified by Chun Ly, 25 November 2017
     - Add optional calib_line keyword option
     - Define dbase and lamp0 variables
     - Use calib_line settings for nsfitcoords and nstransform on telluric dataset
     - Use calib_line settings for nsfitcoords and nstransform on sky-subtracted
       sci dataset
     - Use calib_line settings for nsfitcoords and nstransform on sci OH dataset
    Modified by Chun Ly, 10 December 2017
     - Import glog and call for stdout and ASCII logging
     - glog implementation in prepare, do_flat, do_arcs steps
     - glog implementation in wave_cal, skysub, fitcoords steps
     - glog implementation in combine, extract steps
     - Move up mylogger definition, pass tot0==0 warnings to mylogger
     - Minor bug fix
     - Pass mylogger to normalize_flat()
    Modified by Chun Ly, 18 December 2017
     - Add begin and end QA_clean logging to glog logfile
     - Pass mylogger to iraf_get_subset.check_prefix() calls
     - Pass mylogger to remove_bias_level.run()
     - Pass mylogger to examine_median.run()
    Modified by Chun Ly, 10 January 2018
     - Pass mylogger to wave_cal_script.main()
     - Pass mylogger to iraf_get_subset.main()
    Modified by Chun Ly, 17 January 2018
     - Handle multiple telluric datasets (nsreduce, nsfitcoords, nstransform)
    Modified by Chun Ly, 21 January 2018
     - Handle multiple telluric datasets (nscombine)
    Modified by Chun Ly, 22 January 2018
     - Call mylogger for arc wavelength calibration info
    Modified by Chun Ly, 02 February 2018
     - Bug fix - Fix crash with nsextract with full path given
    Modified by Chun Ly, 20 April 2018
     - Pass mylogger to all examine_median.run() calls
    Modified by Chun Ly, 24 April 2018
     - Handle no telluric data case
     - Handle multiple telluric datasets for remove_bias_level
     - Handle no telluric data case (nsreduce - skysub)
     - Handle no telluric data case (nsfitcoords, nstransform)
     - Handle no telluric data case (nscombine)
     - Handle no telluric data case (nsextract)
    Modified by Chun Ly, 30 April 2018
     - Bug fix: incorrect array name, flatfile_orig -> flatfile
    Modified by Chun Ly, 31 May 2018
     - Call QA_wave_cal.cross_check for arc calibration against OH skylines
    Modified by Chun Ly,  8 June 2018
     - Call QA_wave_cal.residual_wave_cal for both arc/OH with cal = 'arc'
     - Call QA_wave_call.cross_check for OH calibration against arc lines
     - Bug fix: Move residual_wave_cal after cross_check is done
     - Bug fix: OH dataset not ready
    Modified by Chun Ly, 18 June 2018
     - Only log calib_line settings if fitcoords is called
    Modified by Chun Ly, 19 June 2018
     - Call QA_wave_cal.get_database_model
     - Pass function and order to nsfitcoords for arc stack
    Modified by Chun Ly, 20 June 2018
     - Set nsfitcoords xorder fitting
     - Set nsfitcoords x- and y-order fitting in fitcoords
    Modified by Chun Ly, 22 June 2018
     - Call residual_wave_cal for cross check: arc calib for OH skylines
    Modified by Chun Ly, 25 June 2018
     - Check tell_comb files exist before extraction
     - Call iraf_get_subset.check_prefix to check for tfrbnc telluric files
     - Call iraf_get_subset.check_prefix to check for tfrbnc science files
     - Call iraf_get_subset.check_prefix to check for rbnc telluric files
     - Call iraf_get_subset.check_prefix to check for frbnc telluric files
     - Call iraf_get_subset.check_prefix to check for rbnc science files
    Modified by Chun Ly, 26 June 2018
     - Fix mylogger.warn info
     - Bug fix: typo
    Modified by Chun Ly,  6 July 2018
     - Add tell_corr keyword input
     - Call nstelluric
    '''
    
    rawdir = check_path(rawdir) # + on 20/09/2017

    # + on 10/12/2017
    logfile  = rawdir+'reduce.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin run : '+systime())

    # + on 16/05/2017
    tot0 = sum([do_all, prepare, do_flat, do_arcs, wave_cal, skysub,
                fitcoords, combine, extract])
    if tot0 == 0:
        # Mod on 10/12/2017
        mylogger.warn('No GNIRS functions are being called!!!')
        mylogger.warn('Run with do_all=1 or either of these keywords set to 1 : ')
        mylogger.warn('prepare  : Execute nsprepare on all.lis')
        mylogger.warn('do_flat  : Create superflat')
        mylogger.warn('do_arcs  : nsreduce arc data')
        mylogger.warn('wave_cal : Wavelength calibration')
        mylogger.warn('skysub   : Skysubtraction - telluric and science data')
        mylogger.warn('fitcoords: nsfitcoords, nstransform on telluric and science data')
        mylogger.warn('combine  : nscombine telluric and science data')
        mylogger.warn('extract  : nsextract telluric data')
        mylogger.warn('tell_corr : nstelluric on telluric and science data')
        return
    if do_all:
        prepare, do_flat, do_arcs = 1, 1, 1
        wave_cal, skysub = 1, 1
        extract, combine, fitcoords = 1, 1, 1
        tell_corr = 1

    # + on 11/07/2017. Mod on 10/12/2017
    mylogger.info('prepare = %i  do_flat   = %i' % (prepare, do_flat))
    mylogger.info('do_arcs = %i  wave_cal  = %i' % (do_arcs, wave_cal))
    mylogger.info('skysub  = %i  fitcoords = %i' % (skysub, fitcoords))
    mylogger.info('combine = %i  extract   = %i' % (combine, extract))
    if fitcoords:
        mylogger.info('calib_line = %s' % (calib_line)) # + on 25/11/2017
    mylogger.info('tell_corr = %i' % (tell_corr))

    cdir = os.getcwd()+'/' # + on 06/05/2017

    iraf.gemini.nsheaders("gnirs")

    # + on 26/04/2017
    timestamp = systime().replace(':','.')
    logfile   = rawdir+'gnirs_'+timestamp+'.log'
    iraf.gemini.gnirs.logfile = logfile

    mylogger.info("Raw data is located in : %s" % rawdir) # Mod on 10/12/2017

    mylogger.info("GNIRS logfile : "+logfile) # + on 05/05/2017. Mod on 10/12/2017

    # Save reduce.py for each run | + on 05/05/2017
    reduce_file = 'reduce_'+timestamp+'.py'
    mylogger.info("GNIRSLongSlit.reduce script : " + reduce_file) # Mod on 10/12/2017
    os.system('cp -a '+co_filename+' '+rawdir+reduce_file)

    # Check for cleanir files first | Later + on 26/04/2017
    c_files = glob.glob(rawdir+'cN*fits') # Mod on 06/05/2017, 07/06/2017
    if len(c_files) == 0:
        # Mod on 10/12/2017
        mylogger.warn("No cleanir files (cN*fits) available")
        mylogger.warn("Need to execute symlink.run()") # + on 05/05/2017
        mylogger.warn("ABORTING!!!")
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
    OH_obj_list = rawdir+'obj.OH.lis'

    # Handle multiple telluric data is present | + on 17/01/2018
    tell_full_list = glob.glob(rawdir+'telluric?.lis')
    if len(tell_full_list) == 0:
        tell_full_list = [tell_list]

    # + on 17/05/2017
    tell_comb0 = rawdir+'tell_comb.fits'
    obj_comb   = rawdir+'obj_comb.fits'

    # + on 26/04/2017
    all_lis = np.loadtxt(all_list, dtype=type(str)) # Mod on 06/05/2017
    n_all   = len(all_lis)

    # + on 25/11/2017
    if calib_line == 'OH':
        dbase = 'database_OH/'
    if calib_line == 'arc':
        dbase = 'database/'
    lamp0 = 'w'+calib_line+'_stack.fits'

    # Step 1 - Prepare GNIRS data | + on 26/04/2017
    if prepare:
        mylogger.info("Preparing GNIRS data") # Mod on 10/12/2017

        nc_files = glob.glob(rawdir+'ncN*fits') # Mod on 06/05/2017
        n_nc     = len(nc_files)

        if n_nc == n_all:
            mylogger.warn("Files exist! Will not run nsprepare!!") # Mod on 10/12/2017
        else:
            fl_forcewcs = yes
            if n_nc == 0:
                # Mod on 26/04/2017 - Must specify full path
                # Mod on 06/05/2017
                # Mod on 14/09/2017 - Flag CRs
                inimages  = "c@"+all_list
                outimages = rawdir+"nc@"+all_list
                iraf.gnirs.nsprepare(inimages=inimages, rawpath=rawdir,
                                     outprefix='', outimages=outimages,
                                     bpm=bpm, shiftx="INDEF", shifty="INDEF",
                                     fl_forcewcs=fl_forcewcs, fl_correct=no,
                                     fl_cravg=yes)
            else:
                '''
                Warns if files do not exist. Need to implement a way to run
                nsprepare for a subset of data (need to include certain frames)
                '''
                mylogger.warn("The following files do not exist: ") # Mod on 10/12/2017
                iraf_get_subset.main(rawdir, 'nc', all_lis=all_lis, mylogger=mylogger,
                                     silent=True)
                #outfile = 'nc_sub.lis'
                #iraf_get_subset.main(rawdir, 'nc', outfile, all_lis=all_lis)
                #iraf.gnirs.nsprepare(inimages="c@"+outfile, rawpath=rawdir,
                #                     bpm=bpm, shiftx="INDEF", shifty="INDEF",
                #                     fl_forcewcs=fl_forcewcs)
        #endelse

        # + on 14/09/2017
        files0  = np.loadtxt(flat_list, dtype=type(str)).tolist()
        files0 += np.loadtxt(obj_list, dtype=type(str)).tolist()
        # Mod on 24/04/2018
        for _list in tell_full_list:
            if exists(_list):
                files0 += np.loadtxt(_list, dtype=type(str)).tolist()
        files0 = ['nc'+file0 for file0 in files0]

        remove_bias_level.run(rawdir, files0, mylogger=mylogger)
        examine_median.run(rawdir, 'orig', mylogger=mylogger)
        examine_median.run(rawdir, 'bias', mylogger=mylogger)
    #end prepare

    # Step 2 - Create flat | + on 26/04/2017
    if do_flat:
        mylogger.info("Creating super flat") # Mod on 10/12/2017
        flats      = np.loadtxt(flat_list, dtype=type(str)) # Mod on 06/05/2017
        flat_files = ['bnc'+file0 for file0 in flats] # Mod on 06/05/2017 | Mod on 10/05/2017

        tmpflat = rawdir+'tmpflat'
        if not exists(tmpflat):
            if silent == False: mylogger.info('Writing : '+tmpflat) # Mod on 10/12/2017
            np.savetxt(tmpflat, flat_files, fmt='%s')
        else:
            mylogger.warn('File exists!!! : '+tmpflat) # Mod on 10/12/2017

        # tmpflat must not exist prior to call of run() for this call to succeed
        good = computeStatistics(rawdir, flat_files) # Mod on 20/11/2017
        
        if len(good) > 0:
            # Mod on 10/12/2017
            mylogger.info('Flat files to use : ')
            mylogger.info(', '.join(flats[good]))
            np.savetxt(flats_rev, flats[good], fmt='%s')

        # + on 04/05/2017 | Mod on 05/05/2017
        # Mod on 06/05/2017, 18/12/2017
        do_run = iraf_get_subset.check_prefix('rbnc', flats_rev, path=rawdir,
                                              mylogger=mylogger)
        if do_run:
            # Mod on 06/05/2017
            iraf.gnirs.nsreduce(rawdir+'bnc@'+flats_rev, outprefix='',
                                outimages=rawdir+'rbnc@'+flats_rev,
                                fl_sky=no, fl_cut=yes, fl_flat=no, fl_dark=no,
                                fl_nsappwave=no)
        else:
            mylogger.warn("Files exist! Will not run nsreduce!!") # Mod on 10/12/2017

        # + on 05/05/2017
        # Mod on 10/09/2017
        # Mod on 14/09/2017
        if not exists(flatfile):
            iraf.gnirs.nsflat(rawdir+'rbnc@'+flats_rev, flatfile=flatfile)
        else:
            # Mod on 10/12/2017
            mylogger.warn('File exists!!! : '+flatfile)
            mylogger.warn('Will not run nsflat')

        # Generate normalized flat plot | + on 11/07/2017
        # Mod on 10/09/2017, 14/09/2017, 10/12/2017
        normalize_flat(flatfile, rawdir=rawdir, mylogger=mylogger)
    #end do_flat

    # Step 3 : Reduce arcs | + on 05/05/2017
    arcs = np.loadtxt(arc_list, dtype=type(str)) # Mod on 06/05/2017

    if do_arcs:
        mylogger.info("Reducing arc data") # Mod on 10/12/2017
        do_run = iraf_get_subset.check_prefix('rnc', arc_list, path=rawdir,
                                              mylogger=mylogger) # Mod on 18/12/2017
        if do_run:
            # Mod on 06/05/2017
            iraf.gnirs.nsreduce(rawdir+'nc@'+arc_list, outprefix='',
                                outimages=rawdir+'rnc@'+arc_list,
                                fl_sky=no, fl_cut=yes, fl_flat=no,
                                fl_dark=no) #fl_nsappwave=no)
        else:
            # Mod on 10/12/2017
            mylogger.warn('File exists!!!')
            mylogger.warn('Will not run nsreduce on arc data')

        # Combine arcs together | + on 07/11/2017
        arc_stack_file = rawdir+'arc_stack.fits'
        if exists(arc_stack_file):
            # Mod on 10/12/2017
            mylogger.warn('File exists!!!')
            mylogger.warn('Will NOT override '+arc_stack_file)
        else:
            a_files   = [rawdir+'rnc'+file0 for file0 in arcs]
            n_a_files = len(a_files)
            t_hdu     = fits.open(a_files[0]) # + on 08/11/2017
            t_hdr     = t_hdu['SCI'].header
            arcs_arr  = np.zeros( (n_a_files, t_hdr['NAXIS2'],t_hdr['NAXIS1']) )
            arcs_arr2 = np.zeros( (n_a_files, t_hdr['NAXIS2'],t_hdr['NAXIS1']) ) # + on 08/11/2017
            for ii in range(n_a_files):
                arcs_arr[ii]  = fits.getdata(a_files[ii], extname='SCI')
                arcs_arr2[ii] = fits.getdata(a_files[ii], extname='VAR') # + on 08/11/2017
            arcs_avg = np.average(arcs_arr, axis=0)
            # + on 08/11/2017
            arcs_var = np.sqrt(arcs_arr2[0]**2+arcs_arr2[1]**2)/2.0
            t_hdu['SCI'].data = arcs_avg
            t_hdu['VAR'].data = arcs_var
            if silent == False:
                mylogger.info('Writing : '+arc_stack_file) # Mod on 10/12/2017
            t_hdu.writeto(arc_stack_file, output_verify='ignore') # Mod on 08/11/2017

    #end do_arcs

    # Step 4 : Perform wavelength calibration | + on 05/05/2017, Mod on 11/11/2017
    if wave_cal:
        iraf.chdir(rawdir) # + on 20/05/2017
        mylogger.info("Performing interactive wavelength calibration on arc data")

        # + on 12/11/2017
        script_file = 'wave_cal_arc.py'
        if not exists(script_file):
            wave_cal_script.main(rawdir, line_source='arc', mylogger=mylogger)
        else:
            mylogger.info('File exists!!! : '+script_file)
            mylogger.info('Will not override!!!')

        do_run = 0
        if not exists('warc_stack.fits'): do_run = 1
        if do_run:
            mylogger.info("In order to perform interactive calibration, open up")
            mylogger.info("a PyRAF terminal in an anaconda IRAF environment")
            mylogger.info("'cd' into "+rawdir)
            mylogger.info("Execute the following command :")
            mylogger.info("execfile('"+script_file+"')")
            t_out = raw_input("## Hit RETURN when arc wavelength calibration is completed")
        else:
            mylogger.warn('Files exist!!!')
            mylogger.warn('Will not run nswavelength on rnc arc stacked data')

        iraf.chdir(cdir)

        QA_wave_cal.arc_check(rawdir, stack=True) # Mod on 11/11/2017

        # Transform arc data to illustrate expected location of arc lines |  + on 15/06/2017
        iraf.chdir(rawdir)

        # Mod on 12/11/2017
        do_run = 0
        if not exists('farc_stack.fits'): do_run = 1
        if do_run:
            func0, order0 = QA_wave_cal.get_database_model(rawdir, 'arc')
            mylogger.info("Running nsfitcoords on arc stacked data")
            iraf.gnirs.nsfitcoords('arc_stack.fits', outprefix='',
                                   outspectra='farc_stack.fits',
                                   lamp='warc_stack.fits', database='database/',
                                   function=func0, lyorder=order0,
                                   lxorder=QA_wave_cal.xorder)
        else:
            mylogger.warn('Files exist!!!')
            mylogger.warn('Will not run nsfitcoords on rnc arc stacked data')

        # Mod on 12/11/2017
        do_run = 0
        if not exists('tfarc_stack.fits'): do_run = 1
        if do_run:
            mylogger.info("Running nstransform on arc data")
            iraf.gnirs.nstransform('farc_stack.fits', outprefix='',
                                   outspectra='tfarc_stack.fits', database='database/')
        else:
            mylogger.warn('Files exist!!!')
            mylogger.warn('Will not run nstransform on frnc arc stacked data')

        iraf.chdir(cdir)

        # + on 25/05/2017, Mod on 12/11/2017
        QA_wave_cal.arc_check2(rawdir, arcs=arcs, stack=True)

        # Plot arc residuals using arc solution | + on 08/06/2018
        QA_wave_cal.residual_wave_cal(rawdir, dataset='arc', cal='arc')

        # Wavelength calibration with OH skylines
        # + on 13/07/2017
        mylogger.info("Performing non-interactive wavelength calibration on OH data")
        OH_stack.run(rawdir)
        OH_stack.wave_cal(rawdir, cdir)  # + on 17/11/2017
        OH_stack.transform(rawdir) # + on 17/11/2017
        OH_stack.plot_spec(rawdir) # + on 22/11/2017

        # Call cross_check | + on 31/05/2018
        QA_wave_cal.cross_check(rawdir, cdir, 'database/')
        QA_wave_cal.cross_check(rawdir, cdir, 'database_OH/') # + on 08/06/2018

        # Plot OH/arc residuals using arc/OH solution | Mod on 22/06/2018
        QA_wave_cal.residual_wave_cal(rawdir, dataset='arc', cal='OH')
        QA_wave_cal.residual_wave_cal(rawdir, dataset='OH', cal='arc')
    #end wave_cal

    # Step 5a : Sky subtract telluric data | + on 16/05/2017
    if skysub:
        # Mod on 10/09/2017
        fl_flat   = yes
        flatimage = flatfile

        # Mod on 17/01/2018
        for _list in tell_full_list:
            if exists(_list): # Mod on 24/04/2018
                do_run = iraf_get_subset.check_prefix('rbnc', _list, path=rawdir,
                                                      mylogger=mylogger) # Mod on 18/12/2017
                if do_run:
                    mylogger.info("Performing sky subtraction on telluric data, "+_list)
                    iraf.gnirs.nsreduce(rawdir+'bnc@'+_list, outprefix='',
                                        outimages=rawdir+'rbnc@'+_list,
                                        fl_nsappwave=no, fl_sky=yes, fl_flat=fl_flat,
                                        flatimage=flatimage)
                else:
                    mylogger.warn('Files exist!!!')
                    mylogger.warn('Will not run nsreduce on bnc telluric data, '+_list)
            else:
                mylogger.warn('Telluric file does NOT exist : '+_list)

        # Step 5b : Sky subtract science data | + on 16/05/2017
        mylogger.info("Performing sky subtraction on science data")
        do_run = iraf_get_subset.check_prefix('rbnc', obj_list, path=rawdir,
                                              mylogger=mylogger) # Mod on 18/12/2017
        if do_run:
            iraf.gnirs.nsreduce(rawdir+'bnc@'+obj_list, outprefix='',
                                outimages=rawdir+'rbnc@'+obj_list,
                                fl_nsappwave=no, fl_sky=yes,
                                skyimages=rawdir+'bnc@'+sky_list,
                                fl_flat=fl_flat, flatimage=flatimage)
        else:
            mylogger.warn('Files exist!!!')
            mylogger.warn('Will not run nsreduce on bnc sci data')

        examine_median.run(rawdir, 'skysub', mylogger=mylogger) # + on 15/09/2017

        # Run nsreduce without skysubtraction for wavelength check with OH
        # night skylines | + on 01/06/2017
        obj = np.loadtxt(obj_list, dtype=type(str)) # Mod on 06/05/2017
        if not exists(OH_obj_list):
            OH_obj_lists = [file0.replace('.fits','.OH.fits') for file0 in obj]
            if silent == False: mylogger.info('Writing : '+OH_obj_list)
            np.savetxt(OH_obj_list, OH_obj_lists, fmt='%s')
        else:
            mylogger.warn('File exists!!! : '+OH_obj_list)

        do_run = iraf_get_subset.check_prefix('rbnc', OH_obj_list, path=rawdir,
                                              mylogger=mylogger) # Mod on 18/12/2017
        if do_run:
            # Mod on 15/09/2017
            iraf.gnirs.nsreduce(rawdir+'bnc@'+obj_list, outprefix='',
                                outimages=rawdir+'rbnc@'+OH_obj_list,
                                fl_nsappwave=no, fl_sky=no,
                                fl_flat=fl_flat, flatimage=flatimage)

        else:
            mylogger.warn('Files exist!!!')
            mylogger.warn('Will not run nsreduce on bnc sci OH data')

        examine_median.run(rawdir, 'flat', mylogger=mylogger) # + on 15/09/2017
    #end skysub

    # Step 6a : Apply wavelength solution to telluric data | + on 17/05/2017
    if fitcoords:
        # Telluric data
        iraf.chdir(rawdir)

        # Mod on 17/01/2018
        for _list in tell_full_list:
            if exists(_list): # Mod on 24/04/2018
                do_run1 = iraf_get_subset.check_prefix('rbnc', _list,
                                                       path=rawdir, prereq=True)
                if not do_run1:
                    log.warn('rbnc for telluric NOT available!!!')
                    log.warn('Execute reduce.run with skysub=1')
                else:
                    do_run = iraf_get_subset.check_prefix('frbnc', _list, path=rawdir,
                                                          mylogger=mylogger) # Mod on 18/12/2017
                    if do_run:
                        func0, order0 = QA_wave_cal.get_database_model(rawdir, calib_line)
                        mylogger.info("Running nsfitcoords on telluric data, "+_list)
                        iraf.gnirs.nsfitcoords('rbnc@'+_list, outprefix='',
                                               outspectra='frbnc@'+_list,
                                               lamp=lamp0, database=dbase,
                                               function=func0, lyorder=order0,
                                               lxorder=QA_wave_cal.xorder)
                    else:
                        mylogger.warn('Files exist!!!')
                        mylogger.warn('Will not run nsfitcoords on rbnc telluric data, '+_list)

                do_run1 = iraf_get_subset.check_prefix('frbnc', _list,
                                                       path=rawdir, prereq=True)
                if not do_run1:
                    log.warn('frbnc for telluric NOT available!!!')
                else:
                    do_run = iraf_get_subset.check_prefix('tfrbnc', _list, path=rawdir,
                                                          mylogger=mylogger) # Mod on 18/12/2017
                    if do_run:
                        mylogger.info("Running nstransform on telluric data, "+_list)
                        iraf.gnirs.nstransform('frbnc@'+_list, outprefix='',
                                               outspectra='tfrbnc@'+_list,
                                               database=dbase) # Mod on 25/11/2017
                    else:
                        mylogger.warn('Files exist!!!')
                        mylogger.warn('Will not run nstransform on frbnc telluric data, '+_list)
            else:
                mylogger.warn('Telluric file does NOT exist : '+_list)
        #endfor

        # Step 6b : Apply wavelength solution to science data | + on 17/05/2017
        # Science data
        do_run1 = iraf_get_subset.check_prefix('rbnc', obj_list,
                                               path=rawdir, prereq=True)
        if not do_run1:
            log.warn('rbnc for science NOT available!!!')
            log.warn('Execute reduce.run with skysub=1')
        else:
            do_run = iraf_get_subset.check_prefix('frbnc', obj_list, path=rawdir,
                                              mylogger=mylogger) # Mod on 18/12/2017
            if do_run:
                mylogger.info("Running nsfitcoords on science data")
                func0, order0 = QA_wave_cal.get_database_model(rawdir, calib_line)
                iraf.gnirs.nsfitcoords('rbnc@'+obj_list, outprefix='',
                                       outspectra='frbnc@'+obj_list,
                                       lamp=lamp0, database=dbase, function=func0,
                                       lyorder=order0, lxorder=QA_wave_cal.xorder)
            else:
                mylogger.warn('Files exist!!!')
                mylogger.warn('Will not run nsfitcoords on rbnc science data')

            do_run = iraf_get_subset.check_prefix('tfrbnc', obj_list, path=rawdir,
                                                  mylogger=mylogger) # Mod on 18/12/2017
            if do_run:
                mylogger.info("Running nstransform on science data")
                iraf.gnirs.nstransform('frbnc@'+obj_list, outprefix='',
                                    outspectra='tfrbnc@'+obj_list,
                                       database=dbase) # Mod on 25/11/2017
            else:
                mylogger.warn('Files exist!!!')
                mylogger.warn('Will not run nstransform on frbnc science data')

        # Apply wavelength solution to unsubtracted science data
        # + on 02/06/2017
        do_run = iraf_get_subset.check_prefix('frbnc', OH_obj_list, path=rawdir,
                                              mylogger=mylogger) # Mod on 18/12/2017
        if do_run:
            mylogger.info("Running nsfitcoords on science OH data")
            func0, order0 = QA_wave_cal.get_database_model(rawdir, calib_line)
            iraf.gnirs.nsfitcoords('rbnc@'+OH_obj_list, outprefix='',
                                   outspectra='frbnc@'+OH_obj_list,
                                   lamp=lamp0, database=dbase, function=func0,
                                   lyorder=order0, lxorder=QA_wave_cal.xorder)
        else:
            mylogger.warn('Files exist!!!')
            mylogger.warn('Will not run nsfitcoords on rbnc sci OH data')

        # + on 02/06/2017
        do_run = iraf_get_subset.check_prefix('tfrbnc', OH_obj_list, path=rawdir,
                                              mylogger=mylogger) # Mod on 18/12/2017
        if do_run:
            mylogger.info("Running nstransform on science OH data")
            iraf.gnirs.nstransform('frbnc@'+OH_obj_list, outprefix='',
                                   outspectra='tfrbnc@'+OH_obj_list,
                                   database=dbase) # Mod on 25/11/2017
        else:
            mylogger.warn('Files exist!!!')
            mylogger.warn('Will not run nstransform on frbnc sci OH data')

        iraf.chdir(cdir)
        QA_wave_cal.OH_check(rawdir, skysub=False) # + 02/06/2017
        QA_wave_cal.OH_check(rawdir, skysub=True) # + on 25/05/2017
    #end fitcoords

    # Step 7: Combine 2-D spectra | + on 17/05/2017
    if combine:
        # Mod on 21/01/2018
        for tt in range(len(tell_full_list)):
            _list = tell_full_list[tt]
            if exists(_list): # Mod on 24/04/2018
                do_run = iraf_get_subset.check_prefix('tfrbnc', _list,
                                                      path=rawdir, prereq=True)
                if not do_run:
                    log.warn('tfrbnc for telluric NOT available!!!')
                    log.warn('Execute reduce.run with fitcoords=1')
                else:
                    if len(tell_full_list) > 1:
                        tell_comb = tell_comb0.replace('.fits', str(tt)+'.fits')
                    else:
                        tell_comb = tell_comb0
                    if not exists(tell_comb):
                        mylogger.info("Running nscombine on telluric data, "+_list)
                        iraf.gnirs.nscombine(rawdir+'tfrbnc@'+_list, output=tell_comb,
                                             fl_cross=yes, tolerance=0.1)
                    else:
                        mylogger.warn('File exists : '+tell_comb+' !!!')
                        mylogger.warn('Will not run nscombine on tfrbnc telluric data, '+_list)
            else:
                mylogger.warn('Telluric file does NOT exist : '+_list)

        # Mod on 25/06/2018
        do_run = iraf_get_subset.check_prefix('tfrbnc', obj_list,
                                              path=rawdir, prereq=True)
        if not do_run:
            log.warn('tfrbnc for science NOT available!!!')
            log.warn('Execute reduce.run with fitcoords=1')
        else:
            if not exists(obj_comb):
                mylogger.info("Running nscombine on science data")
                iraf.gnirs.nscombine(rawdir+'tfrbnc@'+obj_list, output=obj_comb,
                                     fl_cross=yes, tolerance=0.1)
            else:
                mylogger.warn('File exists : '+obj_comb+' !!!')
                mylogger.warn('Will not run nscombine on tfrbnc science data')

    # Step 8: Extract 1-D spectra | + on 17/05/2017
    if extract:
        # Mod on 01/02/2018, 02/02/2018
        iraf.chdir(rawdir)
        for tt in range(len(tell_full_list)):
            _list = tell_full_list[tt]
            if exists(_list): # Mod on 24/04/2018
                if len(tell_full_list)>1:
                    tell_comb = tell_comb0.replace('.fits', str(tt)+'.fits')
                else:
                    tell_comb = tell_comb0
                if not exists(tell_comb): # + on 25/06/2018
                    log.warn('File does NOT exist : '+tell_comb)
                    log.warn('Execute reduce.run with combine=1 !!!')
                else:
                    outspec = tell_comb.replace('tell','xtell')
                    if not exists(outspec):
                        mylogger.info("Running nsextract on telluric data, "+_list)
                        iraf.gnirs.nsextract(os.path.basename(tell_comb),
                                             outspectra=os.path.basename(outspec),
                                             database='database/')
                    else:
                        mylogger.warn('File exists : '+outspec+' !!!')
                        mylogger.warn('Will not run nsextract on '+_list)
            else:
                mylogger.warn('Telluric file does NOT exist : '+_list)

        iraf.chdir(cdir)

    # Step 9: Determine and apply telluric correction | + on 06/07/2018
    if tell_corr:
        iraf.chdir(rawdir)
        for tt in range(len(tell_full_list)):
            if len(tell_full_list)>1:
                tell_comb = tell_comb0.replace('.fits', str(tt)+'.fits')
            else:
                tell_comb = tell_comb0
            xtell_comb = tell_comb.replace('tell','xtell')

            if not exists(xtell_comb): # + on 25/06/2018
                log.warn('File does NOT exist : '+xtell_comb)
                log.warn('Execute reduce.run with extract=1 !!!')
            else:
                mylogger.info("Running nstelluric with : "+xtell_comb)
                iraf.gnirs.nstelluric(inimages='@nstelluric.lis', cal=xtell_comb,
                                      outspectra='@nstelluric.out.lis',
                                      threshold=0.01, fl_inter=yes,
                                      logfile=logfile)


        iraf.chdir(cdir)

    #os.chdir(cdir) # + on 06/05/2017

    if silent == False: mylogger.info('### End run : '+systime())
#enddef

