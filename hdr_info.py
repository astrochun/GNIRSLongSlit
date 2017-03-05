"""
header_info
===========

Provide description for code here.
"""

import sys, os

from chun_codes import systime

from os.path import exists
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np
import glob

from astropy.table import Table
from astropy import log

def main(path0, silent=False, verbose=True):

    '''
    main() function to obtain information from FITS header and write to
    ASCII file

    Parameters
    ----------
    path0 : str
     Path to FITS file. Must include '/' at the end

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------
    tab0 : astropy.table.Table
     Astropy ASCII table written to [path0]+'info.tbl'
    Notes
    -----
    Created by Chun Ly, 4 March 2017
    '''

    if silent == False: log.info('### Begin main : '+systime())

    fits_files = glob.glob(path0+'N*fits')
    n_files = len(fits_files)

    names0 = ('filename', 'exptime', 'datelabel', 'UT_date', 'obstype',
              'object')
    dtype0 = ('S20', 'f8', 'S30', 'S20', 'S8', 'S100')
    tab0 = Table(names=names0, dtype=dtype0)

    for nn in xrange(n_files):
        basename = os.path.basename(fits_files[nn])
        if verbose == True: log.info('## Reading : '+basename)
        h0 = fits.getheader(fits_files[nn])
        vec0 = [basename, h0['EXPTIME'], h0['DATALAB'],
                h0['DATE-OBS']+'T'+h0['UT'], h0['OBSTYPE'], h0['OBJECT']]
        tab0.add_row(vec0)

    outfile = path0 + 'info.tbl'
    if silent == False: log.info('## Writing : '+outfile)
    asc.write(tab0, outfile, format='fixed_width_two_line')

    if silent == False: log.info('### End main : '+systime())
#enddef
