"""
examine_median
==============

Code intended to check median levels in science data to examine bias
and flatfielding of "large dips"
"""

import sys, os

from chun_codes import systime

from os.path import exists
from astropy.io import fits
from astropy.io import ascii as asc

import numpy as np

import matplotlib.pyplot as plt
import glob

# + on 20/09/2017
from check_path import main as check_path

from astropy import log
from astropy.table import Table

import glog

from astropy.stats import sigma_clipped_stats # + on 29/04/2018

def run(rawdir, style, mylogger=None, silent=False, verbose=True):
    '''
    Main function to generate median plots for various steps in the data
    reduction

    Parameters
    ----------
    rawdir : str
      Path to raw files.

    style : str
      Type of data to plot. Options are:
       'orig': For original data (before any bias removal) 'ncN' files
       'bias': After bias removal 'bncN' files
       'flat': For flattened data 'rbncN' (Will use OH frames since those
                                           are not sky subtracted)
       'skysub': For flattened and skysub data 'rbncN'

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 14 September 2017
    Modified by Chun Ly, 15 September 2017
     - Add if statements for style = 'flat'
     - Add if statements for style = 'skysub'
    Modified by Chun Ly, 20 September 2017
     - Call check_path()
    Modified by Chun Ly, 18 December 2017
     - Implement glog logging, allow mylogger keyword input
    Modified by Chun Ly, 29 April 2018
     - Compute median of background level
     - Remove median from images, plot improved skysubtraction
    Modified by Chun Ly, 30 April 2018
     - Plot aesthetics (axis labels, borders)
     - Write median-remove skysub images
     - Modify plotting to handle style != 'skysub' (single panel)
    Modified by Chun Ly, 5 May 2018
     - Write skysub statistics ASCII table
     - Fix skysub stats ASCII filename
    Modified by Chun Ly, 8 May 2018
     - Force legend in upper right corner
     - Check if median skysub removal was performed (delete rnbc*.v1.fits.gz
       files to rerun everything)
     - Minor fix for different style settings
    '''

    # + on 18/12/2017
    if type(mylogger) == type(None):
        mylog, clog = 0, log
    else:
        mylog, clog = 1, mylogger

    rawdir = check_path(rawdir) # + on 20/09/2017

    if style == 'orig': prefix = 'nc'
    if style == 'bias': prefix = 'bnc'
    if style == 'flat': prefix = 'rbnc' # + on 15/09/2017

    if style == 'skysub': prefix = 'rbnc' # + on 15/09/2017

    if silent == False: clog.info('### Begin run : '+systime())

    infile = rawdir+'obj.lis'
    if silent == False: clog.info('Reading : '+infile)
    files1 = np.loadtxt(infile, dtype=type(str)).tolist()
    files0 = [rawdir+prefix+file1 for file1 in files1]

    # + on 15/09/2017
    if style == 'flat':
        files0 = [file0.replace('.fits','.OH.fits') for file0 in files0]

    if style == 'skysub':
        fig, ax_arr = plt.subplots(nrows=2)
    else:
        fig, ax_arr = plt.subplots()

    if style == 'skysub': # + on 05/05/2018
        arr_mean = np.zeros(len(files0))
        arr_med  = np.zeros(len(files0))
        arr_sig  = np.zeros(len(files0))

    for ii in range(len(files0)):
        # Mod on 08/05/2018
        if style == 'skysub':
            out_file = files0[ii].replace('.fits','.v1.fits.gz')
            if exists(out_file):
                t_file = out_file.replace(rawdir+prefix,'')
                clog.info('Skysub median removed. Reading original skysub file : '+t_file)
                hdu0 = fits.open(out_file) # Mod on 30/04/2018
            else:
                hdu0 = fits.open(files0[ii]) # Mod on 30/04/2018
        else:
            hdu0 = fits.open(files0[ii]) # Mod on 30/04/2018

        im0    = hdu0['SCI'].data
        label0 = files0[ii].replace(rawdir,'')
        med0   = np.median(im0, axis=0)
        if style == 'skysub':
            ax_arr[0].plot(med0, label=label0, linewidth=0.5)
        else:
            ax_arr.plot(med0, label=label0, linewidth=0.5)

        if style == 'skysub': # + on 29/04/2018
            c_mean0, c_med0, c_sig0 = sigma_clipped_stats(med0, sigma=2.0,
                                                          iters=10)
            ax_arr[0].axhline(y=c_med0, color='black', linestyle='--',
                              linewidth=0.5)

            # + on 05/05/2018
            arr_mean[ii] = c_mean0
            arr_med[ii]  = c_med0
            arr_sig[ii]  = c_sig0

            im1  = im0 - c_med0
            med0 = np.median(im1, axis=0)
            ax_arr[1].plot(med0, label=label0, linewidth=0.5)

            # Save gnirs skysub image with v1 suffix | + on 30/04/2018
            clog.info('Saving gnirs skysub file : '+out_file)
            hdu0.writeto(out_file, overwrite=True)

            # Write median-removed skysub image | + on 30/04/2018
            hdu1 = hdu0
            hdu1['SCI'].data = im1
            clog.info('Write median-removed skysub file : '+files0[ii])
            hdu1.writeto(files0[ii], overwrite=True)
        #endif
    #endfor

    if style == 'skysub':
        ax_arr[1].axhline(y=0, color='black', linestyle='--', linewidth=0.5)

        ax_arr[0].legend(fontsize=6, loc='upper right', frameon=False)
        ax_arr[1].legend(fontsize=6, loc='upper right', frameon=False)

        ax_arr[0].xaxis.set_ticklabels([])
        ax_arr[1].set_xlabel('X [pixels]')

        ax_arr[0].set_ylabel('Flux')
        ax_arr[1].set_ylabel('Flux')

        # + on 05/05/2018
        t_arr0 = [files1, arr_mean, arr_med, arr_sig]
        tab0 = Table(t_arr0, names=('File','sky_mean','sky_med','sky_sig'))
        sky_stats_outfile = rawdir+'sky_stats.tbl'
        if silent == False:
            clog.info('Writing : '+sky_stats_outfile)
        asc.write(tab0, sky_stats_outfile, format='fixed_width_two_line')
    else:
        ax_arr.legend(fontsize=6, loc='upper right', frameon=False)
        ax_arr.set_xlabel('X [pixels]')
        ax_arr.set_ylabel('Flux')

    plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99, hspace=0.0)

    out_pdf = rawdir+'median_plot_'+style+'.pdf'
    if silent == False: clog.info('Writing : '+out_pdf)
    fig.savefig(out_pdf)

    if silent == False: clog.info('### End run : '+systime())
#enddef
