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

from astropy import log

def run(rawdir, style, silent=False, verbose=True):
    '''
    Main function to generate median plots for various steps in the data
    reduction

    Parameters
    ----------
    rawdir : str
      Path to raw files. Must end in a '/'

    style : str
      Type of data to plot. Options are:
       'orig': For original data (before any bias removal) 'ncN' files
       'bias': After bias removal 'bncN' files
       'flat': For flattened data 'rbncN' (Will use OH frames since those
                                           are not sky subtracted)

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
    '''

    if style == 'orig': prefix = 'nc'
    if style == 'bias': prefix = 'bnc'
    if style == 'flat': prefix = 'rbnc' # + on 15/09/2017

    if silent == False: log.info('### Begin run : '+systime())

    infile = rawdir+'obj.lis'
    if silent == False: log.info('### Reading : '+infile)    
    files0 = np.loadtxt(infile, dtype=type(str)).tolist()
    files0 = [rawdir+prefix+file0 for file0 in files0]

    # + on 15/09/2017
    if style == 'flat':
        files0 = [file0.replace('.fits','.OH.fits') for file0 in files0]

    fig, ax = plt.subplots()
    
    for ii in range(len(files0)):
        im0    = fits.getdata(files0[ii], extname='SCI')
        label0 = files0[ii].replace(rawdir,'')
        ax.plot(np.median(im0, axis=0), label=label0, linewidth=0.5)

    ax.legend(fontsize=6, frameon=False)

    out_pdf = rawdir+'median_plot_'+style+'.pdf'
    if silent == False: log.info('### Writing : '+out_pdf)
    fig.savefig(out_pdf)

    if silent == False: log.info('### End run : '+systime())
#enddef
