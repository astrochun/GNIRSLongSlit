"""
OH_stack
========

Stack science data to produce 2-D FITS image of OH skylines for wavelength
calibration
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import fits

import numpy as np

import glob

# + on 04/07/2017
from matplotlib import pyplot as plt
from pylab import subplots_adjust
from astropy import log

# + on 04/07/2017
from math import pi
from IQ_plot import gauss1d
from scipy.optimize import curve_fit
#from Zcalbase_gal.observing.locate_em_lines import gaussian_R

# + on 20/09/2017
from check_path import main as check_path

# + on 04/07/2017
from pyraf import iraf #from reduce import iraf

import wave_cal_script # + on 20/11/2017
import get_OH_centers # + on 14/06/2018

import glog # + on 09/01/2018

iraf.gemini(_doprint=0)
iraf.gemini.gnirs(_doprint=0)

log.info("Unlearning tasks")
iraf.gemini.unlearn()
iraf.gemini.gemtools.unlearn()
iraf.gemini.gnirs.unlearn()
iraf.gemini.nsheader('gnirs')

n_sp_pix = 1022 # + on 13/07/2017

co_dirname = os.path.dirname(__file__)
OH_file = co_dirname+'/rousselot2000.dat'

def gaussian(x, mu, sig):
    # + on 04/07/2017
    return 1./(np.sqrt(2.*pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)
#enddef

def gaussian_R(x_arr, lambda_cen, R_spec):
    '''
    Generate an array consisting of a Gaussian profile given the
    spectral resolution

    Parameters
    ----------
    x_arr : array
      An array of wavelengths

    x_lambda : float
      The central wavelength of the Gaussian line

    R_spec : float or double
      Spectral resolution to consider width of emission lines (e.g., R = 3000)

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 4 July 2017
     - Copied from Zcalbase_gal.observing.locate_em_lines
    '''

    t_FWHM = lambda_cen / R_spec # FWHM based on the wavelength of interest
    temp   = gaussian(x_arr, lambda_cen, t_FWHM/(2 * np.sqrt(2*np.log(2))))
    return temp
#enddef

def run(rawdir, silent=False, verbose=False):

    '''
    Main function to combine science data to produce a 2-D FITS image
    containing OH night skylines for wavelength calibration.  Median
    filtering is used to remove detected objects

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
    2-D image containing OH skyline called 'OH_stack.fits'

    Notes
    -----
    Created by Chun Ly, 25 June 2017
    Modified by Chun Ly, 20 September 2017
     - Call check_path()
    Modified by Chun Ly, 16 November 2017
     - Change prefix: rnc to rbnc
    Modified by Chun Ly, 9 January 2018
     - Import glog and call for stdout and ASCII logging
    '''

    # + on 09/01/2018
    logfile  = rawdir+'OH_stack.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin run : '+systime())

    rawdir = check_path(rawdir) # + on 20/09/2017

    obj_list = rawdir + 'obj.lis'
    if silent == False: mylogger.info('Reading : '+obj_list)
    objs = np.loadtxt(obj_list, dtype=type(str))

    rbnc_files = [rawdir+'rbnc'+file0.replace('.fits','.OH.fits') for
                 file0 in objs] # Mod on 16/11/2017

    # Mod on 16/11/2017
    hdu0 = fits.open(rbnc_files[0])
    hdr  = fits.getheader(rbnc_files[0], extname='SCI')
    naxis1 = hdr['NAXIS1']
    naxis2 = hdr['NAXIS2']

    arr0 = np.zeros((len(rbnc_files), naxis2, naxis1)) # Mod on 16/11/2017
    for ii in range(len(rbnc_files)):
        if verbose == True: mylogger.info('Reading : '+obj_list)
        t_data = fits.getdata(rbnc_files[ii], extname='SCI') # Mod on 16/11/2017
        t_med0 = np.median(t_data, axis=0) # Median along columns
        # Remove median along columns
        med_off = np.repeat(t_med0, naxis2).reshape(naxis1,naxis2).transpose()
        arr0[ii] = t_data - med_off

    # Compute median 
    med_arr0 = np.median(arr0, axis=0)

    hdu0['SCI'].data = med_arr0
    hdu0['VAR'].data = np.zeros((naxis2,naxis1))

    outfile = rawdir+'OH_stack.fits'
    stat0   = 'Overwriting' if exists(outfile) else 'Writing'
    if silent == False:
        mylogger.info(stat0+' : '+outfile)
    hdu0.writeto(outfile, output_verify='ignore', overwrite=True)

    if silent == False: mylogger.info('### End run : '+systime())
#enddef

def wave_cal(rawdir, cdir, silent=False, verbose=False):
    '''
    Run gnirs.nswavelength on OH_stack 2-D FITS image for wavelength
    calibration

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
    2-D image with transformation called 'fOH_stack.fits' and 'tfOH_stack.fits'

    Notes
    -----
    Created by Chun Ly, 13 July 2017
    Modified by Chun Ly, 20 September 2017
     - Call check_path()
    Modified by Chun Ly, 16 November 2017
     - Call wave_cal_script to get PyRAF code
    Modified by Chun Ly, 16 November 2017
     - Bug fix: indentation typo with else statement
    Modified by Chun Ly, 20 November 2017
     - Bug fix: Pass in cdir
    Modified by Chun Ly, 9 January 2018
     - Import glog and call for stdout and ASCII logging
    Modified by Chun Ly, 14 June 2018
     - Import and call get_OH_centers
    '''

    # + on 09/01/2018
    logfile  = rawdir+'OH_stack.log'
    mylogger = glog.log0(logfile)._get_logger()

    rawdir = check_path(rawdir) # + on 20/09/2017

    iraf.chdir(rawdir)

    timestamp = systime().replace(':','.')
    logfile   = rawdir+'gnirs_'+timestamp+'.log'
    iraf.gemini.gnirs.logfile = logfile

    get_OH_centers.main(rawdir)

    # + on 16/11/2017
    script_file = 'wave_cal_OH.py'
    if not exists(script_file):
        wave_cal_script.main(rawdir, line_source='OH')
    else:
        # Mod on 09/01/2018
        mylogger.info('File exists!!! : '+script_file)
        mylogger.info('Will not override!!!')

    # + on 16/11/2017
    do_run = 0
    if not exists('wOH_stack.fits'): do_run = 1
    if do_run:
        # Mod on 09/01/2018
        mylogger.info("In order to perform interactive calibration, open up")
        mylogger.info("a PyRAF terminal in an anaconda IRAF environment")
        mylogger.info("'cd' into "+rawdir)
        mylogger.info("Execute the following command :")
        mylogger.info("execfile('"+script_file+"')")
        t_out = raw_input("## Hit RETURN when OH wavelength calibration is completed")
    else:
        # Mod on 09/01/2018
        mylogger.warn('Files exist!!!')
        mylogger.warn('Will not run nswavelength on OH stacked data')

    iraf.chdir(cdir)
#enddef

def transform(rawdir, silent=False, verbose=False):
    '''
    Transform OH_stack 2-D FITS image for wavelength calibration checks

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
    2-D image with transformation called 'fOH_stack.fits' and 'tfOH_stack.fits'

    Notes
    -----
    Created by Chun Ly, 4 July 2017
    Modified by Chun Ly, 20 September 2017
     - Call check_path()
    Modified by Chun Ly, 22 November 2017
     - Add file checking and log.warn calls
    Modified by Chun Ly, 29 November 2017
     - Bug fix: Missing else statement
    Modified by Chun Ly, 9 January 2018
     - Import glog and call for stdout and ASCII logging
    '''

    # + on 09/01/2018
    logfile  = rawdir+'OH_stack.log'
    mylogger = glog.log0(logfile)._get_logger()

    cdir = os.getcwd()+'/' # + on 06/05/2017

    rawdir = check_path(rawdir) # + on 20/09/2017

    iraf.chdir(rawdir)
    mylogger.info("Running nsfitcoords on OH_stack") # Mod on 09/01/2018
    outfile1 = rawdir + 'fOH_stack.fits'
    if not exists(outfile1):
        iraf.gnirs.nsfitcoords('wOH_stack.fits', outprefix='',
                               outspectra='fOH_stack.fits',
                               lamp='wOH_stack.fits', database='database_OH/')
    else:
        # Mod on 09/01/2018
        mylogger.warn('File exists!!! : '+outfile1)
        mylogger.warn('Will not run nsfitcoords on OH stacked data')

    outfile2 = rawdir + 'tfOH_stack.fits'
    if not exists(outfile2):
        iraf.gnirs.nstransform('fOH_stack.fits', outprefix='',
                               outspectra='tfOH_stack.fits',
                               database='database_OH/')
    else:
        # Mod on 09/01/2018
        mylogger.warn('File exists!!! : '+outfile2)
        mylogger.warn('Will not run nstransform on OH stacked data')

    iraf.chdir(cdir)
#enddef


def plot_spec(rawdir, out_pdf='', silent=False, verbose=False):

    '''
    Plot median spectrum of OH skylines

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
    Produces 'OH_spec.pdf'

    Notes
    -----
    Created by Chun Ly, 4 July 2017
     - Overlay Rousselot (2000) model spectrum
    Modified by Chun Ly, 20 September 2017
     - Call check_path()
    Modified by Chun Ly, 9 January 2018
     - Implement glog stdout and ASCII logging
    '''

    # + on 09/01/2018
    logfile  = rawdir+'OH_stack.log'
    mylogger = glog.log0(logfile)._get_logger()

    rawdir = check_path(rawdir) # + on 20/09/2017

    out_pdf = rawdir+'OH_spec.pdf' if out_pdf == '' else rawdir+out_pdf

    im0, hdr = fits.getdata(rawdir+'tfOH_stack.fits', extname='SCI',
                            header=True)

    OH_med0 = np.median(im0, axis=1)

    crval2 = hdr['CRVAL2']
    cd2_2  = hdr['CD2_2']
    crpix  = hdr['CRPIX2']

    x0  = crval2 + cd2_2*np.arange(len(OH_med0))
    x0 /= 1e4 # In microns
    y0  = OH_med0/max(OH_med0)

    fig, ax = plt.subplots()

    if exists(OH_file):
        if silent == False: mylogger.info('Reading : '+OH_file)
        OH_data  = np.loadtxt(OH_file)
        OH_lines = OH_data[:,0] / 1E4
        OH_int   = OH_data[:,1]

        i_max = np.where(OH_med0 == max(OH_med0))[0]
        l_max = x0[i_max[0]]
        p0 = [0.0, max(y0), l_max, 2.0/1e4]
        popt, pcov = curve_fit(gauss1d, x0, y0, p0=p0)
        fwhm   = popt[3]*2*np.sqrt(2*np.log(2))
        R_spec = l_max / fwhm

        OH_spec_mod = np.zeros(len(x0))
        in_rge = np.where((OH_lines >= x0[0]) & (OH_lines <= x0[-1]))[0]
        for ii in range(len(in_rge)):
            idx = in_rge[ii]
            #ax.plot(np.repeat(OH_lines[in_rge[ii]],2), [0,1.1], 'r--')
            temp = OH_int[idx] * gaussian_R(x0, OH_lines[idx], R_spec)
            OH_spec_mod += temp

        ax.plot(x0, y0, 'k-', label='OH_stack')

        OH_spec_mod /= np.max(OH_spec_mod)
        ax.plot(x0, OH_spec_mod, 'r--', alpha=0.75,
                label='OH model (Rousselot 2000)')


    ax.set_xlabel(r'Wavelength ($\mu$m)', fontsize=16)
    ax.set_ylabel('Normalized Flux', fontsize=16)
    ax.set_ylim([0.0, 1.10])
    ax.minorticks_on()
    ax.tick_params(labelsize=14)

    subplots_adjust(left=0.07, right=0.99, bottom=0.07, top=0.99)
    ax.legend(loc='upper right', fontsize=14, frameon=False)
    fig.set_size_inches(11,8)

    fig.savefig(out_pdf)
#enddef

