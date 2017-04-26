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

