"""
header_info
===========

Get information about observations from FITS header
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

# + on 12/03/2017
from . import gnirs_2017a #targets0

import dir_check

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
     Astropy ASCII table written to [path0]+'hdr_info.tbl'

    Notes
    -----
    Created by Chun Ly, 4 March 2017
     - Later re-organized to check for file first
    Modified by Chun Ly, 5 March 2017
     - File exists warning always printed out
     - Include AIRMASS
    Modified by Chun Ly, 23 March 2017
     - Call dir_check.main() to handle multiple date directories
    '''

    if silent == False: log.info('### Begin main : '+systime())

    # + on 23/02/2017
    dir_list = dir_check.main(path0, silent=silent, verbose=verbose)

    if len(dir_list) == 0:
        if silent == False:
            log.info('## No dir found')
        list_path = [path0]
    else:
        if silent == False:
            log.info('## The following date dir found: '+', '.join(dir_list))
        list_path = [path0+a+'/' for a in dir_list]

    for path in list_path:
        outfile = path + 'hdr_info.tbl'

        # Mod later on 04/03/2017 to not overwrite file
        # Mod on 05/03/2017 to always print out this warning
        if exists(outfile):
            log.warning('## File exists : '+outfile)
            log.warning('## Not over-writing!!! ')
        else:
            fits_files = glob.glob(path+'N*fits')
            n_files = len(fits_files)

            # Mod on 05/03/2017 to include airmass
            names0 = ('filename', 'datelabel', 'UT_date', 'obstype', 'object',
                      'exptime', 'airmass', 'grating', 'gratwave', 'filter1',
                      'filter2', 'slit')
            dtype0 = ('S20', 'S30', 'S25', 'S8', 'S100',
                      'f8', 'f8', 'S15', 'f8', 'S10', 'S10', 'S20')
            tab0 = Table(names=names0, dtype=dtype0)

            for nn in xrange(n_files):
                basename = os.path.basename(fits_files[nn])
                if verbose == True: log.info('## Reading : '+basename)
                h0 = fits.getheader(fits_files[nn])
                # Mod on 05/03/2017 to include airmass
                vec0 = [basename, h0['DATALAB'], h0['DATE-OBS']+'T'+h0['UT'],
                        h0['OBSTYPE'], h0['OBJECT'], h0['EXPTIME'],
                        h0['AIRMASS'], h0['GRATING'], h0['GRATWAVE'],
                        h0['FILTER1'], h0['FILTER2'], h0['SLIT']]
                tab0.add_row(vec0)

            if silent == False: log.info('## Writing : '+outfile)
            asc.write(tab0, outfile, format='fixed_width_two_line')
        #endelse
    #endfor

    if silent == False: log.info('### End main : '+systime())
#enddef

def zcalbase_gal_gemini_2017a():
    '''
    Function to run main() on each set of GNIRS 2017A observation set
    to obtain FITS header info

    Parameters
    ----------
    None

    Returns
    -------
    info.tbl : astropy.table.Table
      Astropy ASCII tables outputted to each path

    Notes
    -----
    Created by Chun Ly, 4 March 2017
    Modified by Chun Ly, 12 March 2017
     - global gnirs_2017 use
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = gnirs_2017a # Mod on 12/03/2017
    #targets0 = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        main(path0+target+'/')

#enddef
