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

def main(path0, targets, silent=False, verbose=True):

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
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    Targets, ObsDate, ObsSet = [], [], []
    TotalTime, gratwave, Airmass = [], [], []
    TellStar, TellSet = [], []
    
    for tt in range(len(targets)):
        tt_path = path0+targets[tt]+'/'
        dir_list, list_path = dir_check.main(tt_path, silent=True,
                                             verbose=False)

        for path in list_path:
            Targets.append(targets[tt])

            QA_file = path+'/hdr_info.QA.tbl'
            QA_tab  = asc.read(QA_file, format='fixed_width_two_line')
            idx = [xx for xx in range(len(QA_tab)) if
                   ('obj' in QA_tab['QA'][xx]) or ('sky' in QA_tab['QA'][xx])]

            t_date = QA_tab['UT_date'][idx][0].split('T')[0]
            ObsDate.append(t_date)

            gratwave.append(QA_tab['gratwave'][idx][0])
            exptime = QA_tab['exptime'][idx][0]
            ObsSet.append(str(len(idx))+'x'+str(exptime)+'s')
            TotalTime.append('%.2f' % (len(idx)*exptime/60.0))
            AM0 = QA_tab['airmass'][idx]
            Airmass.append('%.3f-%.3f' % (np.min(AM0),np.max(AM0)))

            t_idx = [xx for xx in range(len(QA_tab)) if
                     ('telluric' in QA_tab['QA'][xx])]
            t_names = list(set(QA_tab['object'][t_idx]))
            print t_names
            t_exptime = QA_tab['exptime'][t_idx][0]
            TellStar.append(QA_tab['object'][t_idx][0])
            TellSet.append(str(len(t_idx))+'x'+str(t_exptime)+'s')

    arr0   = [Targets, ObsDate, ObsSet, TotalTime, gratwave,
              Airmass, TellStar, TellSet]
    names0 = ('Name', 'UT_Date', 'Sequence', 'Int_Time', 'Grating_Wave',
              'Airmass', 'Telluric_Star', 'Tell_Seq')
    tab0 = Table(arr0, names=names0)

    print tab0

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
