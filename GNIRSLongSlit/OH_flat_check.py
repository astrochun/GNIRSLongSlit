"""
OH_flat_check
=============

Code to generate spectrum of OH skyline and examine spatial dependence.
This is to examine if the dip that we see in flats is present in raw data.
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt

from scipy.signal import find_peaks_cwt
from astropy.stats import sigma_clipped_stats

from astropy import log

# + on 20/09/2017
from check_path import main as check_path

def main(rawdir, out_pdf='', silent=False, verbose=True):

    '''
    Main function for OH_flat_check module

    Parameters
    ----------
    rawdir : str
      Path to raw files. Must end in a '/'

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    Produces 'OH_flat_check.pdf'

    Notes
    -----
    Created by Chun Ly, 8 September 2017
    Modified by Chun Ly, 20 September 2017
     - Call check_path()
    '''

    if silent == False: log.info('### Begin main : '+systime())

    rawdir = check_path(rawdir) # + on 20/09/2017

    out_pdf = rawdir+'OH_flat_check.pdf' if out_pdf == '' else rawdir+out_pdf

    infile   = rawdir+'OH_stack.fits'
    if silent == False: log.info('### Reading : '+infile)
    im0, hdr = fits.getdata(infile, extname='SCI', header=True)

    fig, ax = plt.subplots()

    cen = np.int(np.floor(hdr['NAXIS1']/2))
    mid_spec  = np.median(im0[:,cen-5:cen+5], axis=1)
    mid_spec /= np.max(mid_spec)
    med0      = np.median(mid_spec)

    mid_pind0 = np.array(find_peaks_cwt(mid_spec, np.arange(1,3)))

    mid_pind  = mid_pind0[np.where(mid_spec[mid_pind0]-med0 > 0.05)[0]]

    n_peaks = len(mid_pind)

    med_spec = np.zeros((hdr['NAXIS1'], n_peaks))

    for ii in range(n_peaks):
        temp    = im0[mid_pind[ii]-9:mid_pind[ii]+9,:]
        t_spec  = np.median(temp, axis=0)
        mean, median0, std = sigma_clipped_stats(t_spec[250:300],
                                                 sigma=2.0, iters=5)
        t_spec /= median0 #np.max(t_spec)
        med_spec[:,ii] = t_spec
        #plt.plot(t_spec, color='black', linewidth=0.25)

    plt.plot(np.median(med_spec, axis=1), color='black', linewidth=1) #0.25)

    ax.set_xlabel('Spatial X [pixel]', fontsize=16)
    ax.set_ylabel('Normalized Flux', fontsize=16)
    ax.set_ylim([-0.25, 1.5])
    ax.minorticks_on()
    ax.tick_params(labelsize=14)

    fig.set_size_inches(11,8)

    
    if silent == False: log.info('### Writing '+out_pdf)
    fig.savefig(out_pdf, bbox_inches='tight')

    if silent == False: log.info('### End main : '+systime())
#enddef


