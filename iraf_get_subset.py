"""
iraf_get_subset
===============

This code examines existing files and determines which files have not run.
This handles cases where execution crashes or new files are added
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands

import numpy as np

from astropy import log

def main(path, final_prefix, outfile='', all_lis=[], all_file='', 
         silent=False, verbose=True):

    '''
    Check if files exist. If not, create a temporary ASCII file to execute
    IRAF commands for the subset of files that do not exists

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 26 April 2017
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if len(all_lis) == 0 and all_file == '':
        log.warn('## Must specify all_lis or all_file')
        log.warn('## Aborting!!!')
        return
    #endif

    if len(all_lis) == 0:
        all_lis = np.loadtxt(path+'all.lis', dtype=type(str))

    no_files = [file0 for file0 in all_lis if
                exists(path+final_prefix+file0) == False]
    
    if len(no_files) > 0:
        if outfile != '':
            if silent == False: log.info('## Writing : '+path+outfile)
            np.savetxt(path+outfile, no_files, fmt='%s')
        else:
            pref_name = [final_prefix+file0 for file0 in no_files]
            print pref_name

    if silent == False: log.info('### End main : '+systime())
#enddef

def check_prefix(path, final_prefix, input_lis, silent=False, verbose=True):
    '''
    Check if specific files from an input list exist with given prefix

    Parameters
    ----------
    path : str
      Full path to list. Must include '/' at the end

    final_prefix : str
      Files with specific prefix to search for ('rnc', etc).
      This will be added before the filenames in [input_lis]

    input_lis : str
      Filename for input list to check ('arc.lis', 'flat.lis', etc.)

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    do_run : int
      0 - Not all files available
      1 - All files available

    Notes
    -----
    Created by Chun Ly, 5 May 2017
    '''

    if silent == False: log.info('### Begin check_prefix : '+systime())

    input_lis0 = path+input_lis
    if silent == False: log.info('### Reading : '+input_lis0)
    files = np.loadtxt(input_lis0, dtype=type(str))

    f_exist = [file0 for file0 in files if
               exists(path+final_prefix+file0) == True]
    f_noexist = [file0 for file0 in files if
                 exists(path+final_prefix+file0) == False]

    do_run = 0

    if len(f_exist) == 0:
        log.warn('### No files exist')
        do_run = 1
    else:
        if len(f_exist) != len(files):
            log.warn('### Some files do not exist!')
            log.warn(', '.join(f_noexist))
        else:
            log.info('### All files exist!!!')

    if silent == False: log.info('### End check_prefix : '+systime())
    return do_run
#enddef
