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

def main(rawdir, out_pdf='', gz=False, silent=False, verbose=True):

    '''
    Main function to compute and plot statistics

    Parameters
    ----------
    rawdir : str
      Path to raw files.

    gz : boolean
      Indicate using gz-compressed files (mostly for debugging). Default: False

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
     - Plot expected Poissonian level
     - Plot aesthetics (legend, limits, x/ylabel)
     - Compute and plot rms of stacked image
    Modified by Chun Ly, 21 May 2018
     - Bug fix: Mistake in computing RMS for combined 2-D data
     - Write RMS to mylogger; Label 'Stack' RMS on right side of plot
     - Add gz keyword option to use with compressed files
    Modified by Chun Ly, 28 May 2018
     - Add '/' after rawdir if not provided
    '''

    if rawdir[-1] != '/': rawdir = rawdir+'/'

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
            if gz == True: # + on 21/05/2018
                files0 = [t_file+'.gz' for t_file in files0]

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

            #print avg_arr, med_arr, sig_arr

            if out_pdf == '':
                out_pdf = path+'image_stats.pdf'
            else:
                out_pdf = path+out_pdf

            fig, ax_arr = plt.subplots(nrows=2)

            num0 = 1+np.arange(n_files0)
            ax_arr[0].scatter(num0, avg_arr, marker='o', s=50,
                              facecolor='none', edgecolor='k', label='Mean')
            avg_avg = np.average(avg_arr)
            ax_arr[0].axhline(y=avg_avg, c='k', linestyle='dashed')

            ax_arr[0].scatter(num0, med_arr, marker='x', s=25,
                              color='b', label='Median')
            avg_med = np.average(med_arr)
            ax_arr[0].axhline(y=avg_med, c='b', linestyle='dotted')

            ax_arr[0].legend(loc='lower right', fontsize=8)

            ax_arr[1].scatter(num0, sig_arr, marker='o', s=50,
                              facecolor='none', edgecolor='k')
            avg_sig = np.average(sig_arr)
            ax_arr[1].axhline(y=avg_sig, c='k', linestyle='dashed')

            # Expected Poissionan level
            sig_poisson = avg_sig / np.sqrt(n_files0)
            mylogger.info('Expected (Poisson) RMS : %.4f' % sig_poisson) # + on 21/05/2018
            ax_arr[1].axhline(y=sig_poisson, c='b', linestyle='dashed')
            ax_arr[1].text(0, sig_poisson, 'Poisson', color='b', ha='left',
                           va='bottom')

            ax_arr[0].minorticks_on()
            ax_arr[1].minorticks_on()

            ax_arr[0].set_xticklabels([])
            ax_arr[0].set_xlim([-0.25,n_files0+1])
            ax_arr[1].set_xlim([-0.25,n_files0+1])
            ax_arr[1].set_ylabel(r'$\sigma$')
            ax_arr[1].set_xlabel('Frame No.')

            # Plot rms of stacked image
            comb_file = glob.glob(path+'obj_comb.fits')
            if len(comb_file) != 0:
                mylogger.info('Reading : '+comb_file[0])
                comb_data = fits.getdata(comb_file[0], extname='SCI') # Mod on 21/05/2018

                # Mod on 21/05/2018
                c_mean, c_med, c_sig = sigma_clipped_stats(comb_data, sigma=2.0,
                                                           iters=10)
                mylogger.info('Measured RMS in stack : %.4f' % c_sig) # + on 21/05/2018
                ax_arr[1].axhline(y=c_sig, c='g', linestyle='dashed')
                ax_arr[1].text(n_files0+1, c_sig, 'Stack', color='g', ha='right',
                               va='bottom') # Mod on 21/05/2018

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

