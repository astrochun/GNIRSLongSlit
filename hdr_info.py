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

import glog # + on 08/12/2017

def get_offsets(path0, mylogger=None, silent=True, verbose=False):

    '''
    Function to get offsets from FITS header and write it to ASCII file

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
     Astropy ASCII table written to [path0]+'sci_offsets.tbl'

    Notes
    -----
    Created by Chun Ly, 30 May 2017
    Modified by Chun Ly, 10 December 2017
     - Implement glog logging, allow mylogger keyword input
    '''

    # + on 10/12/2017
    if type(mylogger) == type(None):
        mylog, clog = 0, log
    else:
        mylog, clog = 1, mylogger

    if silent == False: log.info('### Begin get_offsets : '+systime())

    dir_list, list_path = dir_check.main(path0, silent=silent, verbose=verbose)

    for path in list_path:
        outfile = path + 'sci_offsets.tbl'

        if exists(outfile):
            # Mod on 10/12/2017
            clog.warning('File exists : '+outfile)
            clog.warning('Not over-writing!!! ')
        else:
            fits_files = np.loadtxt(path+'obj.lis', dtype=type(str))
            fits_files = [path+file0 for file0 in fits_files] # Bug fix
            n_files = len(fits_files)

            names0 = ('filename', 'xoffset', 'yoffset', 'poffset', 'qoffset')
            dtype0 = ('S20', 'f8', 'f8', 'f8', 'f8')
            tab0 = Table(names=names0, dtype=dtype0)

            for nn in xrange(n_files):
                basename = os.path.basename(fits_files[nn])
                if verbose == True: log.info('## Reading : '+basename)
                h0 = fits.getheader(fits_files[nn])
                vec0 = [basename, h0['XOFFSET'], h0['YOFFSET'],
                        h0['POFFSET'], h0['QOFFSET']]
                tab0.add_row(vec0)

            if silent == False: clog.info('Writing : '+outfile)
            asc.write(tab0, outfile, format='fixed_width_two_line')
        #endelse
    #endfor

    if silent == False: clog.info('### End get_offsets : '+systime())
#enddef

def main(path0, silent=False, verbose=True, overwrite=False):

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
    Modified by Chun Ly, 11 May 2017
     - Handle longer filter1 and filter2 FITS values
    Modified by Chun Ly,  8 December 2017
     - Import glog and call for stdout and ASCII logging
    '''

    # Moved up on 10/12/2017
    logfile  = path0+'hdr_info.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin main : '+systime())

    dir_list, list_path = dir_check.main(path0, silent=silent, verbose=verbose)

    # Mod on 23/03/2017
    for path in list_path:
        outfile = path + 'hdr_info.tbl'

        # Mod later on 04/03/2017 to not overwrite file
        # Mod on 05/03/2017 to always print out this warning
        if overwrite == False and exists(outfile):
            # Mod on 08/12/2017
            mylogger.warning('File exists : '+outfile)
            mylogger.warning('Not over-writing!!! ')
        else:
            fits_files = glob.glob(path+'N*fits')
            n_files = len(fits_files)

            # Mod on 05/03/2017 to include airmass
            names0 = ('filename', 'datelabel', 'UT_date', 'obstype', 'object',
                      'exptime', 'airmass', 'grating', 'gratwave', 'filter1',
                      'filter2', 'slit')
            dtype0 = ('S20', 'S30', 'S25', 'S8', 'S100',
                      'f8', 'f8', 'S15', 'f8', 'S20', 'S20', 'S20')
            tab0 = Table(names=names0, dtype=dtype0)

            for nn in xrange(n_files):
                basename = os.path.basename(fits_files[nn])
                if silent == False: mylogger.info('Reading : '+basename)
                h0 = fits.getheader(fits_files[nn])
                # Mod on 05/03/2017 to include airmass
                vec0 = [basename, h0['DATALAB'], h0['DATE-OBS']+'T'+h0['UT'],
                        h0['OBSTYPE'], h0['OBJECT'], h0['EXPTIME'],
                        h0['AIRMASS'], h0['GRATING'], h0['GRATWAVE'],
                        h0['FILTER1'], h0['FILTER2'], h0['SLIT']]
                tab0.add_row(vec0)

            if silent == False: mylogger.info('Writing : '+outfile) # Mod on 08/12/2017
            asc.write(tab0, outfile, format='fixed_width_two_line')
        #endelse
    #endfor

    if silent == False: mylogger.info('### End main : '+systime())
#enddef

def zcalbase_gal_gemini_2017a(offsets=False):
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
    Modified by Chun Ly, 30 May 2017
     - Added option to run get_offsets() with offsets keyword
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = gnirs_2017a # Mod on 12/03/2017
    #targets0 = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        # Mod on 30/05/2017
        if offsets == False:
            main(path0+target+'/')
        else:
            get_offsets(path0+target+'/')
#enddef
