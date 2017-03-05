"""
create_list
===========

Create ASCII .lis for GNIRS reduction pipeline
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
    main() function to sort through data and create individual lists

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
    Created by Chun Ly, 5 March 2017
    '''

    if silent == False: log.info('### Begin main : '+systime())

    infile = path0 + 'hdr_info.tbl'
    if not exists(infile):
        log.warning('### File does not exist : '+infile)
        log.warning('### Exiting!!! '+systime())
        return
    
    if silent == False: log.info('### Reading: '+infile)
    tab0 = asc.read(infile, format='fixed_width_two_line')
    len0 = len(tab0)

    obstype = tab0['obstype']
    object0 = tab0['object']
    filter2 = tab0['filter2']

    i_arc  = [ii for ii in range(len0) if obstype[ii] == 'ARC']
    i_flat = [ii for ii in range(len0) if obstype[ii] == 'FLAT']

    i_tell = [ii for ii in range(len0) if
              (obstype[ii] == 'OBJECT' and
               ('HIP' in object0[ii] or 'HD' in object0[ii]) and
               ('H2_' not in filter2[ii] and 'H_' not in filter2[ii]))]

    i_sci = [ii for ii in range(len0) if
             (obstype[ii] == 'OBJECT' and
              ('HIP' not in object0[ii] and 'HD' not in object0[ii]) and
              ('H2_' not in filter2[ii] and 'H_' not in filter2[ii]))]
    i_sci = np.array(i_sci)

    if len(i_sci) >0:
        i_obj = i_sci[np.arange(0,len(i_sci),2)]
        i_sky = i_sci[np.arange(1,len(i_sci),2)]
        
    if len(i_arc)  == 0: log.warn('## No ARC data found!')
    if len(i_flat) == 0: log.warn('## No FLAT data found!')
    if len(i_tell) == 0: log.warn('## No Telluric data found!')
    if len(i_sci)  == 0: log.warn('## No science data found!')

    prefix = ['arc', 'flat', 'telluric', 'obj', 'sky']
    index  = [i_arc, i_flat, i_tell, i_obj, i_sky]
    zip0   = zip(prefix, index)

    for a,b in zip0:
        outfile = path0+a+'.lis'
        if silent == False: log.info('## Writing : '+outfile)
        np.savetxt(outfile, tab0['filename'][b], fmt='%s')
        #asc.write will not work. Will not produce single column
        #asc.write(tab0[b], outfile, overwrite=True,
        #          format='no_header')
        
    if silent == False: log.info('### End main : '+systime())
#enddef

def zcalbase_gal_gemini_2017a():
    '''
    Function to run main() on each set of GNIRS 2017A observation set
    to obtain ASCII list for GNIRS reduction pipeline

    Parameters
    ----------
    None

    Returns
    -------
    info.tbl : astropy.table.Table
      Astropy ASCII tables outputted to each path

    Notes
    -----
    Created by Chun Ly, 5 March 2017
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        main(path0+target+'/')

#enddef
