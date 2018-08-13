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

from astropy.stats import sigma_clipped_stats
from astropy.convolution import convolve, Box1DKernel
box_kernel = Box1DKernel(5)

from matplotlib.backends.backend_pdf import PdfPages

from . import gnirs_2017a, gnirs_2017b

import dir_check
from IQ_plot import gauss1d
import glog

nrows = 8
ncols = 5

def group(x_cen):
    t_x_cen = x_cen[np.where(x_cen != 0)] # Exclude those without peaks

    ctype  = [''] * x_cen.size #shape[-1]
    mtype  = [''] * x_cen.size #shape[-1]
    labels = [''] * x_cen.size #shape[-1]

    x_cen_s  = np.sort(t_x_cen)
    x_cen_si = np.argsort(t_x_cen)

    ctype0 = ['b','r','g','m','k']
    mtype0 = ['o','s','D','p']

    x_cen_resize = x_cen.reshape(x_cen.size)
    cnt0, cnt1 = 0, 0
    for ii in range(len(x_cen_s)):
        in_rge = np.array([xx for xx in range(x_cen.size) if \
                           (np.absolute(x_cen_resize[xx] - x_cen_s[ii]) <= 3) and
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
            t_val = x_cen_resize[in_rge]
            if len(in_rge) == 1:
                labels[in_rge[0]] = '%.2f' % t_val[0]
            else:
                labels[in_rge[0]] = '%.2f-%.2f' % (min(t_val), max(t_val))

    ctype  = np.array(ctype).reshape(x_cen.shape)
    mtype  = np.array(mtype).reshape(x_cen.shape)
    labels = np.array(labels).reshape(x_cen.shape)
    return ctype, mtype, labels

def group_index(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last: # Part of the group, bump the end
            last = n
        else: # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last # Yield the last group

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
    Modified by Chun Ly, 29 June 2018
     - Write PDF file in each datedir
     - Minor changes to mylogger calls
    Modified by Chun Ly, 31 July 2018
     - Use sigma_clipped_stats; Search for multiple peaks
    Modified by Chun Ly,  1 August 2018
     - Smooth median signal (boxcar); mask for peaks
     - Call group_index(); Require at leat 5 pixels for peak identification
     - Handle multiple peaks
     - Compute number of peaks only once
    Modified by Chun Ly,  2 August 2018
     - Use combine stack for peak identification (more sensitive)
     - Get quick peak centers from combine stack
     - Use peak in med0 if no combine stack or telluric spec
     - Handle peak finding for both telluric and science data
     - Handle multiple peaks in plotting
     - Restrict fitting to within 10 pixels
    Modified by Chun Ly,  3 August 2018
     - Simplify code (x0_bb -> x0)
     - Fix ValueError: Invalid rgba arg ""
     - Fix bugs with n_files and use of 'x' datapoints
     - Use median for x_cen_middle to deal with outliers
     - Compute median using sigma_clipped_stats, exclude outliers from polyfit
     - Switch flag to flag0 to avoid conflict
     - Use peak when combine stack has single line
     - Switch from curve_fit to weighted centering
     - Define p_idx for single peak
    Modified by Chun Ly,  7 August 2018
     - Bug fix: Incorrect x0_diff
     - Handle case when not all skysubtracted frames are used
       (check for obj_rev.lis file)
     - Force integer for index
     - Avoid right edge issue (specific hack for one target)
     - Use DQ array to masked bad pixels and edges
     - Limit weighted computation within expected location
     - Shift x0_max1 for sub-indexing
     - Force x-limit range (zoom-in)
     - Plot aesthetics: ax annotation
    Modified by Chun Ly,  8 August 2018
     - Plot fitting results
     - Bug fix with savefig location, change labeling
    Modified by Chun Ly,  9 August 2018
     - Fix typo with wrong ax subplots use (row determination)
     - Annotate plot with center
     - Plot aesthetics: toplabel, white space, ylim
     - Changes for right center for telluric and single peak science data cases
     - Adjust ax.plot xrange
    Modified by Chun Ly, 10 August 2018
     - Compute median of trace_arr in each bin, best fit polynomial fit
     - Plot best fit to median of trace_arr
     - Update npz savefile with best_fit results
     - Annotate plot with best fit
     - Fix IndexError issue (limit v1 >= 0)
     - Change fig suptitle (just the filename)
    '''

    if rawdir[-1] != '/': rawdir += '/'
    
    logfile  = rawdir+'distort_trace.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin ! ')

    bin_size = 25

    dir_list, list_path = dir_check.main(rawdir, mylogger=mylogger,
                                         silent=silent, verbose=verbose)

    for path in list_path:
        mylogger.info('Working on : '+path.split('/')[-2])
        out_pdf  = path+'distort_trace.pdf'
        out_pdf1 = path+'distort_trace_fit.pdf'
        npz_file = path+'distort_trace.npz'
        obj_file = path+'obj_rev.lis'

        if not exists(npz_file):
            if not exists(obj_file):
                files = glob(path+'tfrbncN????????S????.fits')
            else:
                mylogger.info('File found : '+obj_file)
                npfiles = np.loadtxt(obj_file, dtype=type(str))
                files   = [path+'tfrbnc'+t_file for t_file in npfiles]

                npfilesT = np.loadtxt(path+'telluric.lis', dtype=type(str))
                filesT = [path+'tfrbnc'+t_file for t_file in npfilesT]

                files += filesT
                files.sort()

            n_files = len(files)

            mylogger.info('Number of files found : %i ' % n_files)
            if n_files > 0:
                pdf_pp = PdfPages(out_pdf1)

                hdr0   = fits.getheader(files[0], extname='SCI')
                n_bins = np.int(np.ceil(np.float(hdr0['NAXIS2']))/bin_size)

                x0 = np.arange(hdr0['NAXIS1'])
                y0 = bin_size/2.0 + bin_size * np.arange(n_bins)

                no_c_file = 0
                c_file = glob(path+'obj_comb.fits')
                if len(c_file) == 0:
                    mylogger.warn('Combine frame not found! Using single peak')
                    no_c_file = 1
                    n_peaks = 1
                else:
                    t_c_im = fits.getdata(c_file[0], extname='SCI')
                    c_med0 = np.median(t_c_im, axis=0)

                    # Find peaks | + on 31/07/2018
                    x0_max = np.argmax(c_med0)
                    x0_min = np.argmin(c_med0)

                    idx_mask = np.where((np.absolute(x0-x0_max) >= 10) &
                                        (np.absolute(x0-x0_min) >= 10))[0]
                    sm_c_med0 = convolve(c_med0, box_kernel)

                    t_mean, t_med, t_std = sigma_clipped_stats(sm_c_med0[idx_mask],
                                                               sigma=2, iters=20)
                    idx_det = np.where((sm_c_med0-t_med)/t_std >= 5)[0]

                    list_peak = np.array(list(group_index(idx_det)))
                    peak_idx  = [xx for xx in range(len(list_peak)) if
                                 (list_peak[xx][1]-list_peak[xx][0] >= 5 and
                                  list_peak[xx][0] < 625)] # Avoid right edge issues
                    list_peak = list_peak[peak_idx]

                    n_peaks = len(peak_idx)
                    mylogger.info('Number of peaks found : '+str(n_peaks))

                    peak_ctr = np.zeros(n_peaks)
                    for pp in range(n_peaks):
                        i1, i2 = list_peak[pp][0], list_peak[pp][1]
                        peak_ctr[pp] = i1 + np.argmax(c_med0[i1:i2])

                trace_arr = np.zeros((n_peaks,n_files,n_bins))
                xcen_arr  = np.zeros((n_peaks,n_files))
                fit_arr   = np.zeros((n_peaks,n_files,3))

                for ff in range(n_files):
                    fig, ax = plt.subplots(nrows=nrows, ncols=ncols)

                    t_im0, t_hdr = fits.getdata(files[ff], extname='SCI',
                                               header=True)
                    t_dq = fits.getdata(files[ff], extname='DQ')

                    t_im = np.ma.array(t_im0, mask=t_dq)
                    med0 = np.ma.median(t_im, axis=0)

                    x0_max0 = np.ma.argmax(med0)

                    h_obj = t_hdr['OBJECT']
                    if no_c_file or ('HIP' in h_obj or 'HD' in h_obj):
                        n_peak, use_peak = 1, 1
                    else:
                        n_peak = n_peaks
                        use_peak = 1 if n_peaks == 1 else 0
                        x0_diff = peak_ctr[0] - x0_max0

                    for bb in range(n_bins):
                        row, col = bb / ncols, bb % ncols

                        ty1, ty2 = (0+bb)*bin_size, (1+bb)*bin_size
                        bb_med0 = np.ma.median(t_im[ty1:ty2], axis=0)

                        for pp in range(n_peak):
                            if use_peak == 1:
                                if no_c_file or ('HIP' in h_obj or 'HD' in h_obj):
                                    v1, v2 = x0_max0-15, x0_max0+15
                                else:
                                    v1 = np.max([np.int(list_peak[0][0]-x0_diff-15),0])
                                    v2 = np.int(list_peak[0][1]-x0_diff+15)

                                # print ff, bb, pp, v1, v2, x0_max0

                                x0_max1 = v1+np.ma.argmax(bb_med0[v1:v2])
                                p_idx = np.arange(x0_max1-10,x0_max1+10)
                                # print ff, bb, pp, x0_max0, x0_max1, p_idx[0], p_idx[-1]
                            else:
                                p_idx = np.arange(np.int(list_peak[pp][0]-x0_diff),
                                                  np.int(list_peak[pp][1]-x0_diff))

                            p_med0 = bb_med0[p_idx]

                            x_cen = np.sum(p_med0 * p_idx)/np.sum(p_med0)
                            #p0 = [0.0, y0_max, x0_max, 2.0]
                            #try:
                            #    popt, pcov = curve_fit(gauss1d, x0, bb_med0, p0=p0)
                            #    x_cen = popt[2]
                            #except RuntimeError:
                            #    print 'Runtime error'
                            #    x_cen = p0[2]
                            trace_arr[pp,ff,bb] = x_cen

                            x_off = p_idx[len(p_idx)/2]
                            ax[row,col].plot(p_idx-x_off, p_med0/max(p_med0))
                            ax[row,col].axvline(x=x_cen-x_off)
                            ax[row,col].annotate('%.1f' % x_cen, xy=[0.95,0.95],
                                                 ha='right', va='top', fontsize=9,
                                                 xycoords='axes fraction')
                            if row != nrows-1:
                                ax[row,col].set_xticklabels([])
                            if col != 0:
                                ax[row,col].set_yticklabels([])

                            ax[row,col].set_ylim([-0.05,1.05])
                        #endfor
                    #endfor
                    plt.subplots_adjust(left=0.08, right=0.99, bottom=0.05, top=0.95,
                                        hspace=0.02, wspace=0.02)
                    fig.suptitle(os.path.basename(files[ff]), fontsize=12)
                    fig.savefig(pdf_pp, format='pdf')
                #endfor
                mylogger.info('Writing : '+out_pdf1)
                pdf_pp.close()

                flag0 = np.ones((n_peaks,n_files,n_bins))

                for ff in range(n_files):
                    for pp in range(n_peaks):
                        t_me, t_md, t_s = sigma_clipped_stats(trace_arr[pp,ff],
                                                              sigma=3, iters=10)
                        x_cen_middle      = t_md
                        xcen_arr[pp,ff]   = x_cen_middle
                        trace_arr[pp,ff] -= x_cen_middle

                        diff = trace_arr[pp,ff] - ((y0-512)*0.019)
                        use = np.where(np.absolute(diff) <= 5)[0]
                        if len(use) > 0: flag0[pp,ff,use] = 0

                        fit = np.polyfit(y0[use], trace_arr[pp,ff][use], 2)
                        fit_arr[pp,ff] = fit
                    #endfor
                #endfor

                mylogger.info('Writing : '+npz_file)
                np.savez(npz_file, trace_arr=trace_arr, xcen_arr=xcen_arr,
                         fit_arr=fit_arr, y0=y0, flag0=flag0)
            else:
                mylogger.warn('Files not found !')
        else:
            mylogger.info('Reading : '+npz_file)
            npz = np.load(npz_file)
            trace_arr = npz['trace_arr']
            xcen_arr  = npz['xcen_arr']
            y0        = npz['y0']
            fit_arr   = npz['fit_arr']
            flag0     = npz['flag0']
            n_files   = xcen_arr.shape[1]
            n_peaks   = xcen_arr.shape[0]

        if n_files > 0:
            fig, ax = plt.subplots()
            xlim = [-10,10]

            ctype, mtype, labels = group(xcen_arr)

            x_fit0 = np.zeros(len(y0))
            for bb in range(len(y0)):
                use        = np.where(flag0[:,:,bb] == 0)
                x_fit0[bb] = np.median(trace_arr[use[0],use[1],bb])
            best_fit = np.polyfit(y0, x_fit0, 2)
            mylogger.info('Best fit : [%f, %f, %f]', best_fit[0],
                          best_fit[1], best_fit[2])

            for pp in range(n_peaks):
                for ff in range(n_files):
                    if labels[pp,ff] != '':
                        fc = ctype[pp,ff] if mtype[pp,ff] == 'x' else 'none'
                        ax.scatter(trace_arr[pp,ff,:], y0, marker=mtype[pp,ff],
                                   alpha=0.5, edgecolor=ctype[pp,ff],
                                   facecolor=fc, label=labels[pp,ff])
                        ax.scatter(trace_arr[pp,ff,flag0[pp,ff] == 1],
                                   y0[flag0[pp,ff]==1], marker='x', color='r')

                        pd = np.poly1d(fit_arr[pp,ff])
                        ax.plot(pd(y0), y0, color=ctype[pp,ff], linewidth=0.75,
                                alpha=0.5)

            ax.scatter(x_fit0, y0, marker='o', edgecolor='none',
                       facecolor='black', linewidth=2.0, alpha=0.9)
            best_pd = np.poly1d(best_fit)
            ax.plot(best_pd(y0), y0, color='black', linewidth=2.0, alpha=0.9)
            b_txt = r'x = A y$^2$ + B y + C'+\
                    '\nA = %.3e\nB = %.3e\nC = %.3f' % (best_fit[0],best_fit[1],
                                                        best_fit[2])
            ax.annotate(b_txt, xy=(0.05,0.75), xycoords='axes fraction',
                        ha='left', va='top')

            out_plt   = np.where((trace_arr > xlim[1]) | (trace_arr < xlim[0]))
            n_out_plt = len(out_plt[0])
            flagged   = np.where(flag0 == 1)
            n_flagged = len(flagged[0])

            a_txt = 'N(exclude) = '+str(n_out_plt)+'\nN(flagged) = '+str(n_flagged)
            ax.annotate(a_txt, xy=(0.05,0.95), xycoords='axes fraction', ha='left',
                        va='top')
            ax.legend(loc='lower right')

            ax.set_ylabel('Y [pixels]', fontsize=14)
            ax.set_xlabel('X relative to center [pixels]', fontsize=14)
            ax.set_xlim(xlim) #[-100,100])
            ax.set_ylim([-10,1050])
            fig.suptitle(path)
            ax.minorticks_on()

            fig.set_size_inches(8,8)
            plt.subplots_adjust(left=0.1, right=0.98, bottom=0.07, top=0.98)
            mylogger.info('Writing : '+out_pdf)
            fig.savefig(out_pdf, format='pdf')

            mylogger.info('Updating : '+npz_file)
            np.savez(npz_file, trace_arr=trace_arr, xcen_arr=xcen_arr,
                     fit_arr=fit_arr, y0=y0, flag0=flag0, best_fit=best_fit)
    #endfor

    if silent == False: mylogger.info('### End ! ')
#enddef

def zcalbase_gal_gemini_all():
    '''
    Function to run main() to get spectral distortion corrections for
    ALL GNIRS observations

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 28 June 2018
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = gnirs_2017a + gnirs_2017b
    targets0.sort()

    for target in targets0:
        log.info('## Working on : '+target)
        t_path = path0 + target
        main(t_path)

#enddef

def create_distort_grid(rawdir, silent=False, verbose=True):
    '''
    Create grid (and plot) of distortion for extraction

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 29 June 2018

    Modified by Chun Ly, 10 August 2018
     - Plot best fit

    Modified by Chun Ly, 13 August 2018
     - Rewrite to work for specified rawdir path (instead of all)
    '''

    if rawdir[-1] != '/': rawdir += '/'

    logfile  = rawdir+'distort_trace.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin ! ')

    dir_list, list_path = dir_check.main(rawdir, mylogger=mylogger,
                                         silent=silent, verbose=verbose)

    for path in list_path:
        mylogger.info('Working on : '+path.split('/')[-2])

        fig, ax = plt.subplots(nrows=3)

        npz_files = glob(path+'/????????/distort_trace.npz')
        if len(npz_files) == 0:
            npz_files = glob(path+'distort_trace.npz')

        if len(npz_files) == 0:
            log.warn('No files found!!!')
        else:
            for n_file in npz_files:
                npz = np.load(n_file)
                xcen_arr = npz['xcen_arr']
                fit_arr  = npz['fit_arr']
                best_fit = npz['best_fit']
                n_peaks  = xcen_arr.shape[0]

                for pp in range(n_peaks):
                    ax[0].scatter(xcen_arr[pp], fit_arr[pp,:,2], marker='o',
                                  edgecolor='k', facecolor='none')
                    ax[0].axhline(y=best_fit[2], color='b')
                    ax[1].scatter(xcen_arr[pp], fit_arr[pp,:,1], marker='o',
                                  edgecolor='k', facecolor='none')
                    ax[0].axhline(y=best_fit[1], color='b')
                    ax[2].scatter(xcen_arr[pp], fit_arr[pp,:,0], marker='o',
                                  edgecolor='k', facecolor='none')
                    ax[0].axhline(y=best_fit[0], color='b')

            ax[0].annotate(r'x = A y$^2$ + B y + C', xy=(0.02,0.95),
                           xycoords='axes fraction', ha='left', va='top')
            ax[1].set_ylabel('B')
            ax[2].set_ylabel('A')

            ax[0].set_ylim([-10,0])
            ax[1].set_ylim([-0.05,0.05])
            ax[2].set_ylim([-0.01,0.01])

            ax[0].set_xticklabels([])
            ax[1].set_xticklabels([])

            ax[2].set_xlabel('X [pix]')
            plt.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.1)

        out_pdf = path+'distort_grid.pdf'
        log.info('Writing : '+out_pdf)
        fig.savefig(out_pdf)
#enddef
