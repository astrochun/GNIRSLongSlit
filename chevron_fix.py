"""
chevron_fix
====

Computes median in an 8xN bin and remove bias offsets
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import fits

import numpy as np

import glob

from astropy import log

from scipy.stats import binned_statistic
from astropy.stats import sigma_clipped_stats

def run(files=[''], file_list='', n_rows=1, silent=False, verbose=True):

    '''
    Main function to remove a "chevron" bias pattern noise. This
    noise is diferent from a vertical pattern that can be removed
    by cleanir.py.  The pattern is in individual rows in increments
    of 8 pixels

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 17 June 2017
    '''
    
    if silent == False:
        log.info('### Begin run : '+systime())

    if files[0] == '' and file_list == '':
        log.warn('### No files specified!!!')
        log.warn('Must either specify a list to files or file_list')

    if file_list != '':
        if silent == False: log.info('### Reading : '+file_list)
        files = np.loadtxt(file_list, dtype=type(str)).tolist()

    x0   = xrange(1024)
    bins = 8*np.arange(129)

    for ff in xrange(len(files)):
        hdu0 = fits.open(files[ff])
        im0  = hdu0[1].data

        y_bins = np.int(np.ceil(np.float(im0.shape[0])/n_rows))
        med0_arr = np.zeros(im0.shape)
        for row in range(y_bins):
            mean, median, std = sigma_clipped_stats(im0[row], sigma=2, iters=10, axis=0)
            stat0, edge, num = binned_statistic(x0, im0[row], bins=bins,
                                                statistic='median')

            med0_arr[row] = np.repeat(stat0-median, 8)
        #endfor

        corr_hdu = hdu0
        corr_im0 = im0 - med0_arr
        corr_hdu[1].data = corr_im0

        corr_hdu.writeto(file0.replace('.fits','.bcorr.fits'), overwrite=True)
    #endfor

    if silent == False:
        log.info('### End run : '+systime())
#enddef

