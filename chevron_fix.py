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
    Modified by Chun Ly, 22 June 2017
     - Compute median in each quadrant
     - Apply correction within a certain number of std deviations (handles random noise
       but not regions with low emission from sky)
    Modified by Chun Ly, 15 July 2017
     - Mask pixels above certain sigma threshold
    '''
    
    if silent == False:
        log.info('### Begin run : '+systime())

    if files[0] == '' and file_list == '':
        log.warn('### No files specified!!!')
        log.warn('Must either specify a list to files or file_list')

    if file_list != '':
        if silent == False: log.info('### Reading : '+file_list)
        files = np.loadtxt(file_list, dtype=type(str)).tolist()

    x0   = np.arange(1024)
    bins = 8*np.arange(129)

    for ff in xrange(len(files)):
        hdu0 = fits.open(files[ff])
        im0  = hdu0[1].data
        sh0  = im0.shape
        y_bins = np.int(np.ceil(np.float(im0.shape[0])/n_rows))
        med0_arr = np.zeros(sh0, dtype=np.float64)

        mask0  = np.zeros_like(im0, dtype=np.int) # + on 15/07/2017

        # Compute median for each quadrant | + on 22/06/2017
        median0 = np.zeros(4)
        mean0   = np.zeros(4)
        std0    = np.zeros(4)
        r1 = [sh0[0]/2, sh0[0]/2,        0,        0]
        r2 = [sh0[0],   sh0[0],   sh0[0]/2, sh0[0]/2]
        c1 = [       0, sh0[1]/2,        0, sh0[1]/2]
        c2 = [sh0[1]/2, sh0[1],   sh0[1]/2, sh0[1]]
        for qq in range(4):
            t_im0 = im0[r1[qq]:r2[qq],c1[qq]:c2[qq]]
            mean, median, std = sigma_clipped_stats(t_im0, sigma=2, iters=100)
            mean0[qq]   = mean
            median0[qq] = median
            std0[qq]    = std
            log.info('## q'+str(qq+1)+' %.3f %.3f %3f ' % (mean, median, std))

            # + on 15/07/2017
            src_y, src_x = np.where((t_im0 - median)/std > 2.0)
            mask0[r1[qq]+src_y,c1[qq]+src_x] = 1

        for row in range(y_bins):
            # mean, median, std = sigma_clipped_stats(im0[row], sigma=2, iters=10, axis=0)
            no_mask = np.where(mask0[row] == 0)[0] # + on 15/07/2017

            stat0, edge, num = binned_statistic(x0[no_mask], im0[row][no_mask],
                                                bins=bins, statistic='median')

            # Determine median and std to adopt | + on 22/06/2017
            if row >= sh0[0]/2:
                t_med = np.repeat(median0[0:2], sh0[1]/2)
                t_std = np.repeat(std0[0:2], sh0[1]/2)
            if row < sh0[0]/2:
                t_med = np.repeat(median0[2:4], sh0[1]/2)
                t_std = np.repeat(std0[2:4], sh0[1]/2)

            # Compute offsets | + on 22/06/2017
            med0_arr[row] = np.repeat(stat0, 8) - t_med
            #no_fix = np.where(np.absolute(med0_arr[row]/t_std) >= 5)[0]
            #if len(no_fix) > 0: med0_arr[row,no_fix] = 0.0
        #endfor

        # + on 15/07/2017
        no_fix = np.where(mask0 == 1)
        med0_arr[no_fix] = 0.0

        corr_hdu = hdu0
        corr_im0 = im0 - med0_arr
        corr_hdu[1].data = corr_im0

        corr_hdu.writeto(files[ff].replace('.fits','.bcorr2.fits'), overwrite=True)
        fits.writeto(files[ff].replace('.fits','.fix.fits'), med0_arr, overwrite=True)
        fits.writeto(files[ff].replace('.fits','.mask.fits'), mask0, overwrite=True)

    #endfor

    if silent == False:
        log.info('### End run : '+systime())
#enddef

