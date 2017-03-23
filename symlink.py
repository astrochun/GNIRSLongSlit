"""
symlink
=======

Create symbolic links between raw data for prefix handling of cleanir products

Code checks to see if a cNYYYYMMDDSNNNN.fits exists. If not, it will provide
do a symbolic link
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy.table import Table
from astropy.log import log

def run(path0, silent=False, verbose=True):
    '''
    Main function to run for a given path containing the raw data

    Parameters
    ----------
    path0 : str
     Path to raw FITS file. Must include '/' at the end

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 22 March 2017
    '''
    
    if silent == False: print log.info('### Begin run : '+systime())

    infile0 = path0+'all.lis'
    if silent == False:
        log.info('## Reading : '+infile0)
    files   = np.loadtxt(infile0, dtype=type(str)).tolist()
    n_files = len(files)

    for nn in xrange(n_files):
        c_file = path0+'c'+files[nn]
        if exists(c_file):
            log.info('## File exists : '+c_file)
        else:
            cmd0 = 'ln -fs '+files[nn]+' '+c_file
            print cmd0
    if silent == False: print log.info('### End run : '+systime())
#enddef

