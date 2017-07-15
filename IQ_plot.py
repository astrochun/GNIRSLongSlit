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
from pylab import subplots_adjust

import glob

from astropy.table import Table
from astropy import log

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

# + on 12/03/2017
from . import gnirs_2017a #targets0

import dir_check # + on 11/04/2017

pscale = 0.15 # arcseconds per pixel

# From: https://www.gemini.edu/sciops/telescopes-and-sites/observing-condition-constraints
FWHM_IQ_C = ['20%', '70%', '85%', 'any']
FWHM_IQ_J = [0.40, 0.60, 0.85, 1.40]

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
     - Later modified to shift each median profile and return the array
    Modified by Chun Ly, 11 March 2017
     - Fix interp1d failure when outside bound limits
    Modified by Chun Ly, 11 May 2017
     - Handle RuntimeError with curve_fit
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

    stack0_shift = np.zeros( (len(bins),len(cols)/2) ) # Later + on 10/03/2017

    for bb in xrange(len(bins)):
        temp  = im0_crop[bins[bb]:bins[bb]+bsize-1,:]
        t_med = np.median(temp, axis=0)
        stack0[bb,:] = t_med - np.median(t_med)

        x = np.arange(len(cols))
        y = stack0[bb,:].reshape(len(cols))
        p0 = [0.0, max(y), np.where(y == np.max(y))[0][0], 2.0]

        # Mod on 11/05/2017 to fix RuntimeError bug
        try:
            popt, pcov = curve_fit(gauss1d, x, y, p0=p0)
            fwhm0[bb] = popt[3]*2*np.sqrt(2*np.log(2)) * pscale
        except RuntimeError:
            log.warn('## Optimal parameters not found from curve_fit, %i' % bb)
            fwhm0[bb] = np.nan

        # Later + on 10/03/2017
        x_shift = x - popt[2]
        # Mod on 11/03/2017
        f = interp1d(x_shift, y, bounds_error=False, fill_value=0.0)
        y_shift = f(range(-25,25))
        stack0_shift[bb,:] = y_shift/max(y_shift)

    return bins+bsize/2, fwhm0, stack0_shift

def main(path0='', out_pdf='', check_quality=True, skysub=False, silent=False,
         verbose=True, overwrite=False):
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

    check_quality : boolean
      Check whether data meets IQ requirements. Default: True

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    overwrite : boolean
      Overwrite files if they exists. Default: False

    Returns
    -------
    multi-page PDF plot

    Notes
    -----
    Created by Chun Ly, 10 March 2017
     - Later modified to include check_quality keyword option
     - Later modified to include inset that shows the stacked line profile
    Modified by Chun Ly, 11 April 2017
     - Call dir_check.main() to handle multiple date directories
    Modified by Chun Ly, 13 April 2017
     - Minor bug: Check if file exists first
    Modified by Chun Ly, 10 May 2017
     - Minor bug: When .lis file contains only one entry, problem for
       appending to list
    Modified by Chun Ly, 1 June 2017
     - Added overwrite keyword option to overwrite file. Default is not to
       overwrite .pdf files
     - Bug found: No longer need sky.lis since obj.lis includes all
    Modified by Chun Ly, 6 June 2017
     - Add skysub keyword option to operate on sky-subtracted images
    Modified by Chun Ly, 14 July 2017
     - Fix tick mark locations
     - Fix y limit range for extreme outliers
    '''

    if silent == False: log.info('### Begin main : '+systime())

    # + on 11/04/2017
    dir_list, list_path = dir_check.main(path0, silent=silent, verbose=verbose)

    out_pdf_default = out_pdf

    # Mod on 11/04/2017
    for path in list_path:
        files = []
        file_lis = ['obj.lis','telluric.lis'] # Minor bug fix on 01/06/2017
        for file0 in file_lis:
            # Mod on 13/04/2017
            if exists(path+file0):
                if silent == False: log.info('## Reading : '+path+file0)
                t_files = np.loadtxt(path+file0, dtype=type(str)).tolist()
                # Bug fix - 10/05/2017
                if type(t_files) == str: t_files = [t_files]
                files += t_files
            else:
                if silent == False: log.info('## File not found : '+path+file0)

        files.sort()

        n_files = len(files)

        # Mod on 06/06/2017
        if skysub == True:
            files = ['rnc'+file0 for file0 in files]
            out_pdf = path+'IQ_plot.skysub.pdf' if out_pdf == '' else \
                      path+out_pdf
        else:
            out_pdf = path+'IQ_plot.raw.pdf' if out_pdf == '' else \
                      path+out_pdf

        if overwrite == False and exists(out_pdf):
            log.warn('## File exists!! Will not overwrite '+out_pdf)
        else:
            pp = PdfPages(out_pdf)

            for nn in xrange(n_files):
                if silent == False: log.info('## Reading : '+files[nn])
                hdr0 = fits.getheader(path+files[nn])
                # Mod on 06/06/2017
                if skysub == False:
                    im0 = fits.getdata(path+files[nn])
                else:
                    im0 = fits.getdata(path+files[nn], 'sci')

                bins, fwhm0, stack0_shift = compute_fwhm(im0)

                row = nn % 2
                if row == 0: fig, ax0 = plt.subplots(2, 1)

                good = np.where(fwhm0 >0)[0]
                ax0[row].plot(bins[good], fwhm0[good], marker='o', alpha=0.5,
                              mec='none', mfc='b', linestyle='none', zorder=2)

                ax0[row].get_yaxis().set_tick_params(which='both', right=True,
                                                     width=1, direction='in')
                ax0[row].get_xaxis().set_tick_params(which='both', top=True,
                                                     width=1, direction='in')

                # Compute average line profile from stack0_shift
                # Later + on 10/03/2017
                avg_stack = np.average(stack0_shift, axis=0)
                x0_avg    = np.arange(-25,25)
                if row == 0: axi = fig.add_axes([0.60,0.81,0.25,0.15])
                if row == 1: axi = fig.add_axes([0.60,0.36,0.25,0.15])
                axi.plot(x0_avg*pscale, avg_stack, 'k-')
                axi.set_ylim([-0.4,1.1])
                axi.set_xlim([-2.0,2.0])
                axi.set_xticks(range(-2,2,1))
                axi.set_xlabel('X [arcsec]', fontsize=8)
                axi.tick_params(labelsize=8)
                axi.minorticks_on()
                p0 = [0.0, 1.0, 0.0, 2.0]
                popt, pcov = curve_fit(gauss1d, x0_avg, avg_stack, p0=p0)
                avg_fwhm0 = popt[3]*2*np.sqrt(2*np.log(2)) * pscale
                axi.plot(x0_avg*pscale, gauss1d(x0_avg, *popt), 'r--')
                axi.annotate('FWHM = %.3f"' % avg_fwhm0, [0.50,0.025],
                             xycoords='axes fraction', ha='center', va='bottom')
                ax0[row].axhline(y=avg_fwhm0, linewidth=2, color='r',
                                 linestyle='--', zorder=1)

                # Median FWHM | Later + on 10/03/2017
                med_fwhm0 = np.median(fwhm0[good])
                ax0[row].axhline(y=med_fwhm0, linewidth=2, color='g',
                             linestyle='--', zorder=1)

                # Axes labeling
                ax0[row].set_ylabel('FWHM [arcsec]')
                if row == 1:
                    ax0[row].set_xlabel('Y [PIXELS]')
                else:
                    ax0[row].set_xticklabels([])

                # Annotation
                txt0 = files[nn]+'\nTarget: '+hdr0['OBJECT'] #.split(' ')[0]
                ax0[row].annotate(txt0, [0.025,0.95], xycoords='axes fraction',
                                  ha='left', va='top')

                # Later + on 10/03/2017
                if check_quality:
                    req = hdr0['REQIQ'].replace('-percentile','%')
                    raw = hdr0['RAWIQ'].replace('-percentile','%')

                    i_raw = [ii for ii in range(len(FWHM_IQ_C)) if
                             raw in FWHM_IQ_C[ii]][0]
                    i_req = [ii for ii in range(len(FWHM_IQ_C)) if
                             raw in FWHM_IQ_C[ii]][0]

                    txt0  = 'Req. IQ: %s [%.2f"]\n' % (req, FWHM_IQ_J[i_req])
                    txt0 += 'Raw IQ: %s [%.2f"]\n'  % (raw, FWHM_IQ_J[i_raw])
                    if med_fwhm0 <= FWHM_IQ_J[i_raw]:
                        txt0 += 'PASS'
                    else:
                        if med_fwhm0 <= FWHM_IQ_J[i_raw]*1.25:
                            txt0 += 'USABLE'
                        if med_fwhm0 > FWHM_IQ_J[i_raw]*1.25:
                            txt0 += 'FAIL'

                    ax0[row].annotate(txt0, [0.975,0.05], ha='right',
                                      xycoords='axes fraction', va='bottom')

                # Aesthetics
                ax0[row].set_xlim([0,1050])
                if max(fwhm0[good]) > 3:
                    ax0[row].set_ylim([0,3.0])
                else:
                    ax0[row].set_ylim([min(fwhm0)-0.025,max(fwhm0)+0.075])
                ax0[row].minorticks_on()

                if row == 1 or nn == n_files-1:
                    subplots_adjust(left=0.10, bottom=0.10, top=0.975,
                                    right=0.975, wspace=0.03, hspace=0.05)
                    fig.savefig(pp, format='pdf')
            #endfor
            if silent == False: log.info('## Writing : '+out_pdf)
            pp.close()
        #endelse
        out_pdf = out_pdf_default
    #endfor

    if silent == False: log.info('### End main : '+systime())
#enddef

def zcalbase_gal_gemini_2017a_raw():
    '''
    Function to run main() on each set of GNIRS 2017A observation set
    to obtain PDF plots that illustrate seeing (FWHM) as a function
    of wavelength/pixel

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 10 March 2017
    Modified by Chun Ly, 12 March 2017
     - global gnirs_2017 use
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = gnirs_2017a # Mod on 12/03/2017
    #targets0 = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        main(path0=path0+target+'/', out_pdf='IQ_plot.raw.pdf')
#enddef
