"""
distort_trace
=============

Determine distortion corrections from trace of bright sources
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
from glob import glob

from astropy.table import Table
from astropy import log

from scipy.optimize import curve_fit

from matplotlib.backends.backend_pdf import PdfPages

import dir_check
from IQ_plot import gauss1d
import glog

def group(x_cen):
    ctype  = [''] * len(x_cen)
    mtype  = [''] * len(x_cen)
    labels = [''] * len(x_cen)

    x_cen_s  = np.sort(x_cen)
    x_cen_si = np.argsort(x_cen)

    ctype0 = ['b','r','g','m','k']
    mtype0 = ['o','s','x','D']

    cnt0, cnt1 = 0, 0
    for ii in range(len(x_cen_s)):
        in_rge = np.array([xx for xx in range(len(x_cen)) if \
                           (np.absolute(x_cen[xx] - x_cen_s[ii]) <= 1) and
                           (ctype[xx] == '')])
        #in_rge = np.where((np.absolute(x_cen - x_cen_s[ii]) <= 1) &
        #                  (ctype == ''))[0]
        if len(in_rge) > 0:
            idx = cnt0 % len(ctype0)
            for jj in in_rge:
                ctype[jj] = ctype0[idx]
                mtype[jj] = mtype0[cnt1]
            cnt0 += 1
            if idx == 4: cnt1 += 1
            t_val = x_cen[in_rge]
            if len(in_rge) == 1:
                labels[in_rge[0]] = '%.2f' % t_val[0]
            else:
                labels[in_rge[0]] = '%.2f-%.2f' % (min(t_val), max(t_val))
    return ctype, mtype, labels

def main(rawdir, silent=False, verbose=True):

    '''
    Main function for distort_trace

    Parameters
    ----------
    rawdir : str
      Path to FITS file. Must include '/' at the end

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 27 June 2018
     - Write multi-page PDF file
     - Remove center value for middle of spectra
     - Plot offsets
     - Bug fix for curve_fit (use bb_med0); plot aesthetics (legend)
     - Plot aesthetics (axes labeling), margin adjustments
     - Call group() to get matplotlib markers and colors
     - Write npz file using np.savez and np.load when available
     - Do linear regression fit and save to npz file
    Modified by Chun Ly, 28 June 2018
     - Change linear regression to 2nd order, plot fit
    '''

    if rawdir[-1] != '/': rawdir += '/'
    
    logfile  = rawdir+'distort_trace.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin main ! ')

    bin_size = 25

    dir_list, list_path = dir_check.main(rawdir, mylogger=mylogger,
                                         silent=silent, verbose=verbose)

    out_pdf = rawdir+'distort_trace.pdf'
    pp = PdfPages(out_pdf)
    for path in list_path:
        npz_file = path+'distort_trace.npz'

        if not exists(npz_file):
            files = glob(path+'tfrbncN????????S????.fits')
            n_files = len(files)

            mylogger.info('Number of files found : %i ' % n_files)
            if n_files > 0:
                hdr0      = fits.getheader(files[0], extname='SCI')
                n_bins    = np.int(np.ceil(np.float(hdr0['NAXIS2']))/bin_size)
                trace_arr = np.zeros((n_files,n_bins))
                xcen_arr  = np.zeros(n_files)

                fit_arr   = np.zeros((n_files,3))

                y0 = bin_size/2.0 + bin_size * np.arange(n_bins)

                for ff in range(n_files):
                    t_im = fits.getdata(files[ff], extname='SCI')
                    med0 = np.median(t_im, axis=0)
                    
                    x0_max = np.argmax(med0)
                    for bb in range(n_bins):
                        ty1, ty2 = (0+bb)*bin_size, (1+bb)*bin_size
                        bb_med0 = np.median(t_im[ty1:ty2], axis=0)
                        if bb == 0:
                            x0_bb = np.arange(len(bb_med0))
                        x0_max, y0_max = np.argmax(bb_med0), np.max(bb_med0)
                        p0 = [0.0, y0_max, x0_max, 2.0]
                        try:
                            popt, pcov = curve_fit(gauss1d, x0_bb, bb_med0, p0=p0)
                            x_cen = popt[2]
                        except RuntimeError:
                            print 'Runtime error'
                            x_cen = p0[2]
                        trace_arr[ff,bb] = x_cen
                    #endfor
                    x_cen_middle = trace_arr[ff,n_bins/2]
                    xcen_arr[ff] = x_cen_middle
                    trace_arr[ff] -= x_cen_middle

                    fit = np.polyfit(y0, trace_arr[ff], 2)
                    fit_arr[ff] = fit
                #endfor

                mylogger.info('Writing : '+npz_file)
                np.savez(npz_file, trace_arr=trace_arr, xcen_arr=xcen_arr,
                         fit_arr=fit_arr, y0=y0)
            else:
                mylogger.warn('Files not found !')
        else:
            mylogger.info('Reading : '+npz_file)
            npz = np.load(npz_file)
            trace_arr = npz['trace_arr']
            xcen_arr  = npz['xcen_arr']
            y0        = npz['y0']
            fit_arr   = npz['fit_arr']
            n_files   = len(xcen_arr)

        if n_files > 0:
            fig, ax = plt.subplots()
            xlim = [-10,10]

            ctype, mtype, labels = group(xcen_arr)

            for ff in range(n_files):
                ax.scatter(trace_arr[ff,:], y0, marker=mtype[ff], alpha=0.5,
                           edgecolor=ctype[ff], facecolor='none',
                           label=labels[ff])

                pd = np.poly1d(fit_arr[ff])
                ax.plot(pd(y0), y0, color=ctype[ff], linewidth=0.75, alpha=0.5)

            ax.legend(loc='lower right')

            ax.set_ylabel('Y [pixels]', fontsize=14)
            ax.set_xlabel('X relative to center [pixels]', fontsize=14)
            fig.suptitle(path)
            ax.minorticks_on()

            fig.set_size_inches(8,8)
            plt.subplots_adjust(left=0.1, right=0.98, bottom=0.07, top=0.98)
            fig.savefig(pp, format='pdf')
    #endfor

    pp.close()
    if silent == False: mylogger.info('### End main ! ')
#enddef

