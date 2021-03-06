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

import glog # + on 08/12/2017

def main(path0, silent=False, verbose=True, overwrite=False):

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
    Modified by Chun Ly, 10 April 2017
     - Handle cases when arc, flat, telluric, and sci data are not available
    Modified by Chun Ly, 11 April 2017
     - Handle case when no data for all.lis are available
       (i.e., the observing sequence was canceled)
    Modified by Chun Ly, 16 May 2017
     - Change obj.lis and sky.lis to include both A-B and B-A sets
     - Avoid QA flag of 'sky'
     - Fix minor bug
    Modified by Chun Ly, 17 May 2017
     - Minor fix for i_sky to handle skipping of frames for i_obj
    Modified by Chun Ly, 31 May 2017
     - Minor fix for when no science data is available
    Modified by Chun Ly, 1 June 2017
     - Added overwrite keyword option to overwrite file
    Modified by Chun Ly, 8 December 2017
     - Import glog and call for stdout and ASCII logging
    Modified by Chun Ly, 10 December 2017
     - Adjust glog logging to a default file in path0, move up mylogger
       definition
     - Minor mylogger text changes
     - Pass mylogger to dir_check.main()
    Modified by Chun Ly, 12 January 2018
     - Handle multiple telluric datasets
    Modified by Chun Ly, 15 January 2018
     - Handle multiple telluric datasets (cont'd) - Get indices for each dataset
    Modified by Chun Ly, 16 January 2018
     - Handle multiple telluric datasets (cont'd) - Update prefix and index lists
    Modified by Chun Ly, 17 January 2018
     - Create telluric.lis even for multi-telluric case
    Modified by Chun Ly, 23 April 2018
     - Handle no telluric data case
    '''

    logfile  = path0+'create_list.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin main : '+systime())

    # + on 23/03/2017
    dir_list, list_path = dir_check.main(path0, mylogger=mylogger,
                                         silent=silent, verbose=verbose)

    # Mod on 23/03/2017
    for path in list_path:
        infile = path + 'hdr_info.tbl'
        if not exists(infile):
            mylogger.warning('File does not exist : '+infile)
            mylogger.warning('Exiting!!! '+systime())
            return
    
        if silent == False: mylogger.info('Reading: '+infile)
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

        # Identify when multiple sets of telluric data is available
        # Broken based on filename | + on 12/01/2018, Mod on 15/01/2018
        # tel_name = list(set(np.array(object0[i_tell])))
        multi_tell = 0

        # Mod on 23/04/2018
        if len(i_tell) != 0:
            seq_tell = np.array([np.int(str0.replace('.fits','')[-4:]) for
                                 str0 in tab0['filename'][i_tell]])
            seq_diff = seq_tell[1:] - seq_tell[0:-1]
            # +1 to factor first skip for seq_diff
            seq_break = [ii+1 for ii in range(len(seq_diff)) if seq_diff[ii] != 1]
            if len(seq_break) > 0:
                multi_tell = 1
                i_tell0 = []
                mylogger.warn('Multiple telluric star detected.')
                mylogger.warn('Will split into separate files.')
                ss_start = [0] + seq_break # relative to i_tell
                ss_end   = seq_break + [len(i_tell)]
                for ss in range(len(ss_start)):
                    ss_idx = i_tell[ss_start[ss]:ss_end[ss]]
                    i_tell0.append(ss_idx)
            else:
                mylogger.info('Only one telluric star detected.')

        i_sci = [ii for ii in r0 if
                 (obstype[ii] == 'OBJECT' and
                  ('HIP' not in object0[ii] and 'HD' not in object0[ii]) and
                  ('H2_' not in filter2[ii] and 'H_' not in filter2[ii]))]
        i_sci = np.array(i_sci)

        # Note: Assumes that dithering pattern is ABA'B'
        # Mod on 16/05/2017
        # Mod on 17/05/2017 - Minor bug if things skip for sci frames
        # Mod on 31/05/2017 - Minor bug if i_sci is empty
        if len(i_sci) > 0:
            i_obj = i_sci
            i_off = [1, -1] * (len(i_sci)/2)
            if len(i_sci) % 2 == 1: i_off.append(-1) # Odd number correction
            i_sky = i_sci[np.arange(len(i_sci))+np.array(i_off)]
            # i_sky = [a+b for a,b in zip(i_sci,i_off)]

        # Mod on 10/04/2017
        prefix, index = [], []
        if len(i_arc)  == 0:
            mylogger.warn('No ARC data found!')
        else:
            prefix.append('arc')
            index.append(i_arc)

        if len(i_flat) == 0:
            mylogger.warn('No FLAT data found!')
        else:
            prefix.append('flat')
            index.append(i_flat)

        if len(i_tell) == 0:
            mylogger.warn('No Telluric data found!')
        else:
            # Mod on 16/01/2018, 17/01/2018
            # Always create telluric.lis for all tellurics
            prefix.append('telluric')
            index.append(i_tell)
            if multi_tell:
                for ii in range(len(i_tell0)):
                    prefix.append('telluric'+str(ii+1))
                    index.append(i_tell0[ii])

        if len(i_sci)  == 0:
            mylogger.warn('No science data found!')
        else:
            prefix.append('obj')
            prefix.append('sky')
            index.append(i_obj)
            index.append(i_sky)

        zip0   = zip(prefix, index)

        QA = ['N/A'] * len0 # Later + on 05/03/2017

        for a,b in zip0:
            if a != 'sky':
                for idx in b: QA[idx] = a # Mod on 16/05/2017
            outfile = path+a+'.lis'
            if overwrite == False and exists(outfile):
                mylogger.warn('File exists! Will not write '+outfile+'!!')
            else:
                if silent == False: mylogger.info('Writing : '+outfile)
                np.savetxt(outfile, tab0['filename'][b], fmt='%s')
                #asc.write will not work. Will not produce single column
                #asc.write(tab0[b], outfile, overwrite=True,
                #          format='no_header')

        # Later + on 05/03/2017 | Mod on 11/04/2017
        i_all    = [ii for ii in r0 if QA[ii] != 'N/A']
        if len(i_all) > 0:
            outfile0 = path+'all.lis'
            if overwrite == False and exists(outfile0):
                mylogger.warn('File exists! Will not write all.lis!!')
            else:
                if silent == False: mylogger.info('Writing : '+outfile0)
                np.savetxt(outfile0, tab0['filename'][i_all], fmt='%s')
        else: mylogger.warn('Will not write all.lis!!')

        # Later + on 05/03/2017
        col0 = Column(QA, name='QA')
        tab0.add_column(col0)

        # Later + on 05/03/2017
        outfile2 = infile.replace('.tbl', '.QA.tbl')
        if silent == False:
            if overwrite == False and exists(outfile2):
                mylogger.warn('File exists! Will not write '+outfile2+'!!')
            else:
                if not exists(outfile2):
                    mylogger.info('Writing : '+outfile2)
                else: mylogger.info('Overwriting : '+outfile2)
                asc.write(tab0, outfile2, format='fixed_width_two_line',
                          overwrite=True)

    if silent == False: mylogger.info('### End main : '+systime())
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
