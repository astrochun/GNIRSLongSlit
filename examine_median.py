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

import numpy as np

import matplotlib.pyplot as plt
import glob

# + on 20/09/2017
from check_path import main as check_path

from astropy import log

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
    files0 = np.loadtxt(infile, dtype=type(str)).tolist()
    files0 = [rawdir+prefix+file0 for file0 in files0]

    # + on 15/09/2017
    if style == 'flat':
        files0 = [file0.replace('.fits','.OH.fits') for file0 in files0]

    fig, ax = plt.subplots()
    
    for ii in range(len(files0)):
        im0    = fits.getdata(files0[ii], extname='SCI')
        label0 = files0[ii].replace(rawdir,'')
        med0   = np.median(im0, axis=0)
        ax.plot(med0, label=label0, linewidth=0.5)

        if style == 'skysub': # + on 29/04/2018
            c_mean0, c_med0, c_sig0 = sigma_clipped_stats(med0, sigma=2.0,
                                                          iters=10)
            ax.axhline(y=c_med0, color='black', linestyle='--', linewidth=0.5)

    ax.legend(fontsize=6, frameon=False)

    out_pdf = rawdir+'median_plot_'+style+'.pdf'
    if silent == False: clog.info('Writing : '+out_pdf)
    fig.savefig(out_pdf)

    if silent == False: clog.info('### End run : '+systime())
#enddef
