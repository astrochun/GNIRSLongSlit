"""
IQ_plot
=======

Generate plots that uses bright alignment star to determine IQ
"""

import sys, os

from chun_codes import systime

from os.path import exists
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import glob

from astropy.table import Table
from astropy import log

from scipy.optimize import curve_fit

import astropy.units as u
pscale = 0.15 #Arcsec per pixel

from pylab import subplots_adjust

def gauss1d(x, a0, a, x0, sigma):
    return a0 + a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def compute_fwhm(im0):
    '''
    Determines FWHM as a function of y-pixel or wavelength. Uses the brightest
    source in the longslit spectra, and computes a median profile with a bin
    size of 50 pixels and fits the median profile with a Gaussian using
    curve_fit()

    Parameters
    ----------
    im0 : numpy array
      Array from fits.getdata() of the FITS image


    Returns
    -------
    bins : numpy array
     Array of y pixels

    fwhm0 : numpy array
     Array of FWHM values in each bin

    Notes
    -----
    Created by Chun Ly, 10 March 2017
    '''

    # First find the peak of the nearby star
    med_arr = np.median(im0, axis=0) # Average along columns
    i_max   = np.where(med_arr == np.max(med_arr))[0][0]
    cols    = range(i_max-50,i_max+50)
    
    im0_crop = im0[:,cols]
    naxis2   = im0.shape[0]

    bsize = 50
    bins  = np.arange(0, naxis2, bsize)

    stack0 = np.zeros( (len(bins),len(cols)) )
    fwhm0 = np.zeros(len(bins))

    for bb in xrange(len(bins)):
        temp = im0_crop[bins[bb]:bins[bb]+bsize-1,:]
        stack0[bb,:] = np.median(temp, axis=0)

        x = np.arange(len(cols))
        y = stack0[bb,:].reshape(len(cols))
        p0 = [0.0, max(y), np.where(y == np.max(y))[0][0], 2.0]
        popt, pcov = curve_fit(gauss1d, x, y, p0=p0)
        fwhm0[bb] = popt[3]*2*np.sqrt(2*np.log(2)) * pscale

    return bins+bsize/2, fwhm0

def main(path0='', out_pdf='', silent=False, verbose=True):
    '''
    main() function to compute natural seeing (image quality) from bright
    alignment star

    Parameters
    ----------
    path0 : str
      Full path to where output PDF and FITS file are located. Must end
      with a '/'

    out_pdf : str
      Filename for output PDF. Do NOT include full path

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    multi-page PDF plot

    Notes
    -----
    Created by Chun Ly, 10 March 2017
    '''

    if silent == False: log.info('### Begin main : '+systime())

    files = []
    file_lis = ['obj.lis','sky.lis','telluric.lis']
    for file0 in file_lis:
        if silent == False: log.info('## Reading : '+path0+file0)    
        t_files = np.loadtxt(path0+file0, dtype=type(str)).tolist()
        files += t_files
    files.sort()
    
    n_files = len(files)

    out_pdf = path0+'IQ_plot.pdf' if out_pdf == '' else path0+out_pdf

    pp = PdfPages(out_pdf)

    for nn in xrange(n_files):
        if silent == False: log.info('## Reading : '+files[nn])    
        hdr0 = fits.getheader(path0+files[nn])
        im0  = fits.getdata(path0+files[nn])

        bins, fwhm0 = compute_fwhm(im0)

        row = nn % 2
        if row == 0: fig, ax_arr = plt.subplots(2, 1)

        good = np.where(fwhm0 >0)[0]
        ax_arr[row].plot(bins[good], fwhm0[good], marker='o', alpha=0.5,
                         mec='none', mfc='b', linestyle='none')
        ax_arr[row].set_ylabel('FWHM [arcsec]')
        if row == 1:
            ax_arr[row].set_xlabel('Y [PIXELS]')
        else:
            ax_arr[row].set_xticklabels([])

        txt0 = files[nn]+'\n'+hdr0['OBJECT']
        ax_arr[row].annotate(txt0, [0.025,0.95], xycoords='axes fraction',
                             ha='left', va='top')
        
        ax_arr[row].set_xlim([0,1050])
        ax_arr[row].set_ylim([min(fwhm0)-0.05,max(fwhm0)+0.05])
        ax_arr[row].minorticks_on()
        if row == 1 or nn == n_files-1:
            subplots_adjust(left=0.10, bottom=0.10, top=0.975, right=0.975,
                            wspace=0.03, hspace=0.05)
            fig.savefig(pp, format='pdf')
    #endfor

    if silent == False: log.info('## Reading : '+out_pdf)
    pp.close()

    if silent == False: log.info('### End main : '+systime())
#enddef

def zcalbase_gal_gemini_2017a_raw():
    '''
    Function to run main() on each set of GNIRS 2017A observation set
    to obtain PDF plots of each *raw* FITS image for visual examination

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 6 March 2017
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        main(path0=path0+target+'/', out_pdf='IQ_plot.raw.pdf')
#enddef
