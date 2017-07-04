"""
OH_stack
========

Stack science data to produce 2-D FITS image of OH skylines for wavelength
calibration
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import fits

import numpy as np

import glob

# + on 04/07/2017
from matplotlib import pyplot as plt
from pylab import subplots_adjust

from astropy import log

# + on 04/07/2017
from pyraf import iraf #from reduce import iraf
iraf.gemini(_doprint=0)
iraf.gemini.gnirs(_doprint=0)

log.info("Unlearning tasks")
iraf.gemini.unlearn()
iraf.gemini.gemtools.unlearn()
iraf.gemini.gnirs.unlearn()
iraf.gemini.nsheader('gnirs')

def run(rawdir, silent=False, verbose=False):

    '''
    Main function to combine science data to produce a 2-D FITS image
    containing OH night skylines for wavelength calibration.  Median
    filtering is used to remove detected objects

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
    2-D image containing OH skyline called 'OH_stack.fits'

    Notes
    -----
    Created by Chun Ly, 25 June 2017
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    obj_list = rawdir + 'obj.lis'
    if silent == False: log.info('### Reading : '+obj_list)
    objs = np.loadtxt(obj_list, dtype=type(str))

    rnc_files = [rawdir+'rnc'+file0.replace('.fits','.OH.fits') for
                 file0 in objs]

    hdu0 = fits.open(rnc_files[0])
    hdr  = fits.getheader(rnc_files[0], extname='SCI')
    naxis1 = hdr['NAXIS1']
    naxis2 = hdr['NAXIS2']

    arr0 = np.zeros((len(rnc_files), naxis2, naxis1))
    for ii in range(len(rnc_files)):
        if verbose == True: log.info('### Reading : '+obj_list)
        t_data = fits.getdata(rnc_files[ii], extname='SCI')
        t_med0 = np.median(t_data, axis=0) # Median along columns
        # Remove median along columns
        med_off = np.repeat(t_med0, naxis2).reshape(naxis1,naxis2).transpose()
        arr0[ii] = t_data - med_off

    # Compute median 
    med_arr0 = np.median(arr0, axis=0)

    hdu0['SCI'].data = med_arr0
    hdu0['VAR'].data = np.zeros((naxis2,naxis1))

    outfile = rawdir+'OH_stack.fits'
    stat0   = 'Overwriting' if exists(outfile) else 'Writing'
    if silent == False:
        log.info('## '+stat0+' : '+outfile)
    hdu0.writeto(outfile, output_verify='ignore', overwrite=True)

    if silent == False: log.info('### End run : '+systime())
#enddef

def transform(rawdir, silent=False, verbose=False):
    '''
    Transform OH_stack 2-D FITS image for wavelength calibration checks

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
    2-D image containing OH skyline called 'OH_stack.fits'

    Notes
    -----
    Created by Chun Ly, 4 July 2017
    '''

    cdir = os.getcwd()+'/' # + on 06/05/2017

    iraf.chdir(rawdir)
    log.info("## Running nsfitcoords on OH_stack")
    iraf.gnirs.nsfitcoords('wOH_stack.fits', outprefix='',
                           outspectra='fOH_stack.fits',
                           lamp='wOH_stack.fits', database='database_OH/')

    iraf.gnirs.nstransform('fOH_stack.fits', outprefix='',
                           outspectra='tfOH_stack.fits',
                           database='database_OH/')
    iraf.chdir(cdir)
#enddef

