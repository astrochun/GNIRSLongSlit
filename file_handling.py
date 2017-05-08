"""
file_handling
=============

A set of functions to copy and remove files

"""

import sys, os

from chun_codes import systime
from chun_codes import match_nosort
from chun_codes import chun_crossmatch

from os.path import exists
import commands
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
import glob

from astropy.table import Table
from astropy import log

def cp_files(outdir, path, final_prefix, input_lis, silent=False, verbose=True):

    '''
    Copy files into outdir for certain IRAF functions

    Parameters
    ----------
    outdir : str
      Full path for where files are temporarily stored.
      Must include '/' at the end

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

    Notes
    -----
    Created by Chun Ly, 7 May 2017
    '''
    
    if silent == False: log.info('### Begin cp_files : '+systime())
    
    if silent == False: log.info('### Reading : '+input_lis)
    files = np.loadtxt(input_lis, dtype=type(str))

    files0 = [path+final_prefix+file0 for file0 in files]

    cmd0 = 'cp -a '+' '.join(files0)+' '+outdir
    if silent == False: log.info('### '+cmd0)
    os.system(cmd0)

    if silent == False: log.info('### End cp_files : '+systime())
#enddef

def rm_files(path, final_prefix, input_lis, silent=False, verbose=True):

    '''
    Delete files in the current directory

    Parameters
    ----------
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

    Notes
    -----
    Created by Chun Ly, 7 May 2017
    '''
    
    if silent == False: log.info('### Begin rm_files : '+systime())
    
    input_lis0 = path+input_lis
    if silent == False: log.info('### Reading : '+input_lis0)
    files = np.loadtxt(input_lis0, dtype=type(str))

    files0 = [final_prefix+file0 for file0 in files]

    cmd0 = 'rm '+' '.join(files0)
    if silent == False: log.info('### '+cmd0)
    os.system(cmd0)

    if silent == False: log.info('### End rm_files : '+systime())
#enddef
