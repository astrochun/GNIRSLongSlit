"""
cleanir_script
====

Provide description for code here.
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

import cleanir

co     = __file__
co_dir = os.path.dirname(co)+'/'
cmd0   = co_dir+'cleanir.py'
        
def run(path0, clean_file='', out_script='', silent=False, verbose=True):

    '''
    Run cleanir.py for a set of files

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
    Created by Chun Ly, 7 March 2017
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    if clean_file == '': clean_file = 'clean.lis'
    
    if silent == False: log.info('## Reading : '+path0+clean_file)
    files = np.loadtxt(path0+clean_file, dtype=type(str)).tolist()

    if out_script == '': out_script = path0+'run_cleanir.sh'

    if silent == False: log.info('## Writing : '+out_script)
    f = open(out_script, 'w')
    for ii in xrange(len(files)):
        cmd1 = cmd0+' -afqo '+path0+'c'+files[ii]+' '+path0+files[ii]
        f.write(cmd1+'\n')

    f.close()

    if silent == False: log.info('### End run : '+systime())
#enddef

def zcalbase_gal_gemini_2017a():
    '''
    Function to run run() on each set of GNIRS 2017A observation set
    to generate a script to run on the command line

    Parameters
    ----------
    None

    Returns
    -------
    info.tbl : astropy.table.Table
      Astropy ASCII tables outputted to each path

    Notes
    -----
    Created by Chun Ly, 7 March 2017
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = ['DEEP05', 'DEEP06', 'DEEP07'] #, 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        run(path0+target+'/', out_script=path0+'run_cleanir.'+target+'.sh')

#enddef
