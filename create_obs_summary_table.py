"""
create_obs_summary_table
====

Provide description for code here.
"""

import sys, os

from chun_codes import systime

from os.path import exists
from astropy.io import ascii as asc

import numpy as np

import glob

from . import gnirs_2017a
import dir_check

from astropy.table import Table
from astropy import log

def get_telluric_info(QA_tab, t_idx, t_names):
    ## + on 25/04/2017

    if len(t_idx) == 1:
        telstar   = QA_tab['object'][t_idx][0]
        t_exptime = QA_tab['exptime'][t_idx][0]
        telset    = str(len(t_idx))+'x'+str(t_exptime)+'s'

        AM0   = QA_tab['airmass'][t_idx]
        telAM = '%.3f-%.3f' % (np.min(AM0),np.max(AM0))
    else:
        telset = []
        telAM  = []
        for nn in range(len(t_names)):
            idx2      = [xx for xx in range(len(QA_tab)) if
                         (t_names[nn] in QA_tab['object'][xx])]
            nidx      = list(set(t_idx) & set(idx2))
            t_exptime = QA_tab['exptime'][nidx][0]
            telset.append(str(len(nidx))+'x'+str(t_exptime)+'s')

            AM0   = QA_tab['airmass'][nidx]
            telAM.append('%.3f-%.3f' % (np.min(AM0),np.max(AM0)))

        telstar = ','.join(t_names)
        telset  = ','.join(telset)
        telAM   = ','.join(telAM)
    return telstar, telset, telAM
#enddef

def main(path0, targets, outfile=None, silent=False, verbose=True):

    '''
    Generate ASCII file summarizing observations

    Parameters
    ----------
    path0 : str
      Parent path for all files

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 25 April 2017
    Modified by Chun Ly, 5 May 2017
     - Handle overwriting file
    Modified by Chun Ly, 3 June 2017
     - Bug fix: Check if hdr_info.QA.tbl exists
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    Targets, ObsDate, ObsSet = [], [], []
    TotalTime, gratwave, Airmass = [], [], []
    TellStar, TellSet, TellAM = [], [], [] # Later Mod on 25/04/2017

    for tt in range(len(targets)):
        tt_path = path0+targets[tt]+'/'
        dir_list, list_path = dir_check.main(tt_path, silent=True,
                                             verbose=False)

        cnt = 0
        for path in list_path:
            # Mod on 25/04/2017
            txt = targets[tt] if cnt == 0 else '...'
            Targets.append(txt) #targets[tt])

            QA_file = path+'/hdr_info.QA.tbl'
            # Mod on 03/06/2017
            if not exists(QA_file):
                log.warn('## File not found! '+QA_file)

                ObsDate.append('N/A')
                gratwave.append('N/A')

                ObsSet.append('N/A')
                TotalTime.append('N/A')
                Airmass.append('N/A')

                TellStar.append('N/A')
                TellSet.append('N/A')
                TellAM.append('N/A')
            else:
                QA_tab  = asc.read(QA_file, format='fixed_width_two_line')

                # All science targets
                idx = [xx for xx in range(len(QA_tab)) if
                       ('obj' in QA_tab['QA'][xx]) or
                       ('sky' in QA_tab['QA'][xx])]

                # Later Mod on 25/04/2017
                tab_ref = QA_tab[idx][0]
                t_date  = tab_ref['UT_date'].split('T')[0]
                exptime = tab_ref['exptime']

                ObsDate.append(t_date)
                gratwave.append(tab_ref['gratwave'])

                ObsSet.append(str(len(idx))+'x'+str(exptime)+'s')
                TotalTime.append('%.2f' % (len(idx)*exptime/60.0))

                AM0 = QA_tab['airmass'][idx]
                Airmass.append('%.3f-%.3f' % (np.min(AM0),np.max(AM0)))

                t_idx = [xx for xx in range(len(QA_tab)) if
                         ('telluric' in QA_tab['QA'][xx])]
                t_names = list(set(QA_tab['object'][t_idx]))

                telstar, telset, telAM = get_telluric_info(QA_tab, t_idx, t_names)
                TellStar.append(telstar)
                TellSet.append(telset)
                TellAM.append(telAM) # Later + on 25/04/2017
                cnt += 1
            #endelse
        #endfor
    #endfor

    arr0   = [Targets, ObsDate, ObsSet, TotalTime, gratwave,
              Airmass, TellStar, TellSet, TellAM]
    names0 = ('Name', 'UT_Date', 'Sequence', 'Int_Time', 'Grating_Wave',
              'Airmass', 'Telluric_Star', 'Telluric_Seq', 'Telluric_AM')
    tab0 = Table(arr0, names=names0)

    print tab0

    if outfile == None: outfile = path0+'obs_summary.txt'

    # Mod on 06/05/2017
    if silent == False:
        stat0 = 'Overwriting : ' if exists(outfile) else 'Writing : '
        log.info(stat0+outfile)
    asc.write(tab0, output=outfile, format='fixed_width_two_line', overwrite=True)

    if silent == False: log.info('### End main : '+systime())
#enddef

def zcalbase_gal_gemini_2017a():
    '''
    Function to run main() to get summary table for GNIRS 2017A observations

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

    targets0 = gnirs_2017a

    main(path0, targets0)

#enddef
