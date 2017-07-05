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

# + on 04/07/2017
from pyraf import iraf #from reduce import iraf
iraf.gemini(_doprint=0)
iraf.gemini.gnirs(_doprint=0)

log.info("Unlearning tasks")
iraf.gemini.unlearn()
iraf.gemini.gemtools.unlearn()
iraf.gemini.gnirs.unlearn()
iraf.gemini.nsheader('gnirs')

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
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    obj_list = rawdir + 'obj.lis'
    if silent == False: log.info('### Reading : '+obj_list)
    objs = np.loadtxt(obj_list, dtype=type(str))

    rnc_files = [rawdir+'rnc'+file0.replace('.fits','.OH.fits') for
                 file0 in objs]

    hdu0 = fits.open(rnc_files[0])
    hdr  = fits.getheader(rnc_files[0], extname='SCI')
    naxis1 = hdr['NAXIS1']
    naxis2 = hdr['NAXIS2']

    arr0 = np.zeros((len(rnc_files), naxis2, naxis1))
    for ii in range(len(rnc_files)):
        if verbose == True: log.info('### Reading : '+obj_list)
        t_data = fits.getdata(rnc_files[ii], extname='SCI')
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
        log.info('## '+stat0+' : '+outfile)
    hdu0.writeto(outfile, output_verify='ignore', overwrite=True)

    if silent == False: log.info('### End run : '+systime())
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
    '''

    cdir = os.getcwd()+'/' # + on 06/05/2017

    iraf.chdir(rawdir)
    log.info("## Running nsfitcoords on OH_stack")
    iraf.gnirs.nsfitcoords('wOH_stack.fits', outprefix='',
                           outspectra='fOH_stack.fits',
                           lamp='wOH_stack.fits', database='database_OH/')

    iraf.gnirs.nstransform('fOH_stack.fits', outprefix='',
                           outspectra='tfOH_stack.fits',
                           database='database_OH/')
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
    '''

    out_pdf = rawdir+'OH_spec.pdf' if out_pdf == '' else rawdir+out_pdf

    im0, hdr = fits.getdata(rawdir+'tfOH_stack.fits', extname='SCI',
                            header=True)

    OH_med0 = np.median(im0, axis=1)

    crval2 = hdr['CRVAL2']
    cd2_2  = hdr['CD2_2']
    crpix  = hdr['CRPIX2']

    x0  = crval2 + cd2_2*np.arange(len(OH_med0))
    x0 /= 1e4 # In microns
    fig, ax = plt.subplots()
    ax.plot(x0, OH_med0/max(OH_med0), 'k-', label='OH_stack')

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

