"""
remove_bias_level
====

Remove bias levels for each quadrant of GNIRS detector
"""

import sys, os

from chun_codes import systime

from os.path import exists
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np
import glob

# + on 20/09/2017
from check_path import main as check_path

from astropy.table import Table
from astropy import log

import glog

def run(rawdir, files0=[''], mylogger=None, silent=False, verbose=False):

    '''
    Main function to execute removal of bias in individual GNIRS images

    Parameters
    ----------
    rawdir : str
      Path to raw files.

    files0 : list. Optional
      Name of files. No need to specify full path as rawdir is used

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 14 September 2017
     - Change files0 to not include full path
    Modified by Chun Ly, 20 September 2017
     - Call check_path()
    Modified by Chun Ly, 18 December 2017
     - Implement glog logging, allow mylogger keyword input
     - Fix indentation for End info log
    Modified by Chun Ly, 23 April 2018
     - Write tab_outfile only once (outside for loop)
    '''

    # + on 18/12/2017
    if type(mylogger) == type(None):
        mylog, clog = 0, log
    else:
        mylog, clog = 1, mylogger

    if silent == False: clog.info('### Begin run : '+systime())

    rawdir = check_path(rawdir) # + on 20/09/2017

    if files0[0] == '':
        files0 = glob.glob(rawdir+'ncN????????S????.fits')
        
    n_files0 = len(files0)

    bias_offset0 = np.zeros(n_files0)
    bias_offset1 = np.zeros(n_files0)
    bias_offset2 = np.zeros(n_files0)
    bias_offset3 = np.zeros(n_files0)

    for nn in range(n_files0):
        outfile = rawdir+files0[nn].replace('ncN','bncN') # Mod on 14/09/2017
        if exists(outfile):
            clog.warn('Will not overwrite file : '+outfile)
        else:
            if verbose == True: clog.info('Reading : '+outfile)
            hdu = fits.open(rawdir+files0[nn]) # Mod on 14/09/2017

            if ('BIAS_FIX' in hdu['SCI'].header) == False:
                naxis1 = hdu['SCI'].header['NAXIS1']
                naxis2 = hdu['SCI'].header['NAXIS2']
                qxsize = naxis1 / 2
                qysize = naxis2 / 2
                
                im0 = hdu['SCI'].data
                
                q1x1,q1x2, q1y1,q1y2 =   0,160,    qysize,naxis2
                q2x1,q2x2, q2y1,q2y2 = 864,naxis1, qysize,naxis2
                q3x1,q3x2, q3y1,q3y2 =   0,160,         0,qysize
                q4x1,q4x2, q4y1,q4y2 = 864,naxis1,      0,qysize

                bias_UL = np.median(im0[q1y1:q1y2,q1x1:q1x2])
                bias_UR = np.median(im0[q2y1:q2y2,q2x1:q2x2])
                bias_LL = np.median(im0[q3y1:q3y2,q3x1:q3x2])
                bias_LR = np.median(im0[q4y1:q4y2,q4x1:q4x2])

                bias_offset0[nn] = bias_UL
                bias_offset1[nn] = bias_UR
                bias_offset2[nn] = bias_LL
                bias_offset3[nn] = bias_LR

                im0[qysize:naxis2,     0:qxsize] -= bias_UL
                im0[qysize:naxis2,qxsize:naxis1] -= bias_UR
                im0[     0:qysize,     0:qxsize] -= bias_LL
                im0[     0:qysize,qxsize:naxis1] -= bias_LR

                if verbose == True: clog.info('Writing : '+outfile)
                hdu['SCI'].data = im0
                hdu['SCI'].header['BIAS_FIX'] = 'Yes'
                hdu.writeto(outfile, overwrite=True, output_verify='ignore')
            #endif
        #endelse
    #endfor

    # Mod on 23/04/2018
    arr0   = [files0, bias_offset0, bias_offset1, bias_offset2, bias_offset3]
    names0 = names=('filename','bias_UL','bias_UR','bias_LL','bias_LR')
    tab0   = Table(arr0, names=names0)

    # Mod on 23/04/2018
    tab_outfile = rawdir+'bias_levels.tbl'
    if silent == False: clog.info('Writing : '+tab_outfile)
    asc.write(tab0, tab_outfile, format='fixed_width_two_line')

    if silent == False: clog.info('### End run : '+systime())
#enddef

