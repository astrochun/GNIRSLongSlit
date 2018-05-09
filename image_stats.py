"""
image_stats
===========

Plot statistics of different 2-D sky-subtracted images and compare against
stacked results
"""

import sys, os

from os.path import exists

from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy.table import Table
from astropy import log

import dir_check
import glog

from astropy.stats import sigma_clipped_stats

def main(rawdir, out_pdf='', silent=False, verbose=True):

    '''
    Main function to compute and plot statistics

    Parameters
    ----------
    rawdir : str
      Path to raw files.

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 8 May 2018
     - Change output PDF filename
     - Plot averages of statistic measurements
    '''

    logfile  = rawdir+'image_stats.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin run ! ')

    dir_list, list_path = dir_check.main(rawdir, mylogger=mylogger,
                                         silent=silent, verbose=verbose)

    out_pdf_default = out_pdf

    for path in list_path:
        file_lis = path+'obj.lis'
        if not exists(file_lis):
            if silent == False: mylogger.info('File not found : '+file_lis)
        else:
            if silent == False: mylogger.info('Reading : '+file_lis)
            t_files = np.loadtxt(file_lis, dtype=type(str)).tolist()
            files0  = [path+'tfrbnc'+t_file for t_file in t_files]

            n_files0 = len(files0)
            avg_arr  = np.zeros(n_files0)
            med_arr  = np.zeros(n_files0)
            sig_arr  = np.zeros(n_files0)
            
            for nn in range(n_files0):
                im0 = fits.getdata(files0[nn], extname='SCI')

                c_mean, c_med, c_sig = sigma_clipped_stats(im0, sigma=2.0,
                                                           iters=10)
                avg_arr[nn] = c_mean
                med_arr[nn] = c_med
                sig_arr[nn] = c_sig
            #endfor

            print avg_arr, med_arr, sig_arr

            if out_pdf == '':
                out_pdf = path+'image_stats.pdf'
            else:
                out_pdf = path+out_pdf

            fig, ax_arr = plt.subplots(nrows=2)

            num0 = np.arange(n_files0)
            ax_arr[0].scatter(num0, avg_arr, marker='o', s=50,
                              facecolor='none', edgecolor='k')
            avg_avg = np.average(avg_arr)
            ax_arr[0].axhline(y=avg_avg, c='k', linestyle='dashed')

            ax_arr[0].scatter(num0, med_arr, marker='x', s=25,
                              color='b')
            avg_med = np.average(med_arr)
            ax_arr[0].axhline(y=avg_med, c='b', linestyle='dotted')

            ax_arr[1].scatter(num0, sig_arr, marker='o', s=50,
                              facecolor='none', edgecolor='k')
            avg_sig = np.average(sig_arr)
            ax_arr[1].axhline(y=avg_sig, c='k', linestyle='dashed')
            
            #fig, ax_arr = plt.subplots(ncols=2)
            #ax_arr[0].hist(avg_arr, bins=5, align='mid', color='b',
            #               linestyle='solid', alpha=0.5, edgecolor='b',
            #               histtype='step', label='Average')
            #ax_arr[0].hist(med_arr, bins=5, align='mid', # color='g',
            #               linestyle='dashed', edgecolor='k',
            #               histtype='step', label='Median')
            #ax_arr[1].hist(sig_arr, bins=5, align='mid', color='k',
            #               linestyle='solid', alpha=0.5, edgecolor='k',
            #               histtype='step')

            mylogger.info('Writing : '+out_pdf)
            fig.set_size_inches(8,4)
            fig.savefig(out_pdf, bbox_inches='tight')
            plt.close()

        out_pdf = out_pdf_default
    #endfor

    
    if silent == False: mylogger.info('### End run ! ')
#enddef

