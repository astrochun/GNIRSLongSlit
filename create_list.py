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

from astropy.table import Table, Column
from astropy import log

# + on 12/03/2017
from . import gnirs_2017a #targets0

import dir_check

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
    ASCII files containing the filenames for each set
     obj.lis      - List of science data frames
     sky.lis      - List of science data frames to be used for sky-subtraction
     arc.lis      - List of Arclamp data frames
     flat.lis     - List of Flat data frames
     telluric.lis - List of Telluric data frames
     all.lis      - List of all obj, sky, arc, flat, and telluric data frames

    tab0 : astropy.table.Table
     Astropy ASCII table with QA flag written to [path0]+'hdr_info.QA.tbl'

    Notes
    -----
    Created by Chun Ly, 5 March 2017
     - Later modified to output ASCII table (hdr_info.QA.tbl) containing a
       column called 'QA' with information of how each FITS file was
       classified. This is to enable quick check that each dataset if
       properly classified
     - Later modified to include all.lis output
     - Later modified to define r0 (for slight speed improvement
    Modified by Chun Ly, 23 March 2017
     - Call dir_check.main() to handle multiple date directories
    '''

    if silent == False: log.info('### Begin main : '+systime())

    # + on 23/03/2017
    dir_list = dir_check.main(path0, silent=silent, verbose=verbose)

    # + on 23/03/2017
    if len(dir_list) == 0:
        if silent == False:
            log.info('## No dir found')
        list_path = [path0]
    else:
        if silent == False:
            log.info('## The following date dir found: '+', '.join(dir_list))
        list_path = [path0+a+'/' for a in dir_list]

    # Mod on 23/03/2017
    for path in list_path:
        infile = path + 'hdr_info.tbl'
        if not exists(infile):
            log.warning('### File does not exist : '+infile)
            log.warning('### Exiting!!! '+systime())
            return
    
        if silent == False: log.info('### Reading: '+infile)
        tab0 = asc.read(infile, format='fixed_width_two_line')
        len0 = len(tab0)
        r0   = xrange(len0)

        obstype = tab0['obstype']
        object0 = tab0['object']
        filter2 = tab0['filter2']

        i_arc  = [ii for ii in r0 if obstype[ii] == 'ARC']
        i_flat = [ii for ii in r0 if obstype[ii] == 'FLAT']

        i_tell = [ii for ii in r0 if
                  (obstype[ii] == 'OBJECT' and
                   ('HIP' in object0[ii] or 'HD' in object0[ii]) and
                   ('H2_' not in filter2[ii] and 'H_' not in filter2[ii]))]

        i_sci = [ii for ii in r0 if
                 (obstype[ii] == 'OBJECT' and
                  ('HIP' not in object0[ii] and 'HD' not in object0[ii]) and
                  ('H2_' not in filter2[ii] and 'H_' not in filter2[ii]))]
        i_sci = np.array(i_sci)

        # Note: Assumes that dithering pattern is ABA'B'
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

        QA = ['N/A'] * len0 # Later + on 05/03/2017

        for a,b in zip0:
            for idx in b: QA[idx] = a
            outfile = path+a+'.lis'
            if silent == False: log.info('## Writing : '+outfile)
            np.savetxt(outfile, tab0['filename'][b], fmt='%s')
            #asc.write will not work. Will not produce single column
            #asc.write(tab0[b], outfile, overwrite=True,
            #          format='no_header')

        # Later + on 05/03/2017
        i_all    = [ii for ii in r0 if QA[ii] != 'N/A']
        outfile0 = path+'all.lis'
        if silent == False: log.info('## Writing : '+outfile0)
        np.savetxt(outfile0, tab0['filename'][i_all], fmt='%s')

        # Later + on 05/03/2017
        col0 = Column(QA, name='QA')
        tab0.add_column(col0)

        # Later + on 05/03/2017
        outfile2 = infile.replace('.tbl', '.QA.tbl')
        if silent == False:
            if not exists(outfile2):
                log.info('## Writing : '+outfile2)
            else: log.info('## Overwriting : '+outfile2)
        asc.write(tab0, outfile2, format='fixed_width_two_line', overwrite=True)

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
    Modified by Chun Ly, 12 March 2017
     - global gnirs_2017 use
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = gnirs_2017a # Mod on 12/03/2017
    #targets0 = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        main(path0+target+'/')

#enddef
