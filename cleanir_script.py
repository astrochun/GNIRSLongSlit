"""
cleanir_script
====

Creates UNIX-based .sh script to run cleanir.py
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

# + on 12/03/2017
from . import gnirs_2017a #targets0
import dir_check # + on 15/05/2017

co     = __file__
co_dir = os.path.dirname(co)+'/'
cmd0   = co_dir+'cleanir.py'
        
def run(path0, clean_file='', out_script='', silent=False, verbose=True):

    '''
    Create a .sh script to run cleanir.py for a set of files

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
    Modified by Chun Ly, 15 May 2017
     - Call dir_check.main() to handle multiple date directories
    '''
    
    if silent == False: log.info('### Begin run : '+systime())

    if clean_file == '': clean_file = 'clean.lis'

    # + on 15/05/2017
    dir_list, list_path = dir_check.main(path0, silent=silent, verbose=verbose)

    # Mod on 15/05/2017
    for date,path in zip(dir_list,list_path):
        clean_file0 = path+clean_file
        if not exists(clean_file0):
            log.warn('### File does not exist!!!')
            log.info('### '+clean_file0)
        else:
            if silent == False: log.info('## Reading : '+clean_file0)
            files = np.loadtxt(clean_file0, dtype=type(str)).tolist()

            out_script0 = path+'run_cleanir.sh' if out_script == '' \
                          else out_script

            if len(dir_list) != 0:
                out_script0 = out_script0.replace('.sh', '.'+date+'.sh')

            if silent == False:
                stat0 = 'Overwriting' if exists(out_script0) else 'Writing'
                log.info('## '+stat0+' : '+out_script0)
            f = open(out_script0, 'w')
            for ii in xrange(len(files)):
                cmd1 = cmd0+' -afqo '+path+'c'+files[ii]+' '+path+files[ii]
                f.write(cmd1+'\n')

            f.close()
        #endelse
    #endfor

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
    Modified by Chun Ly, 12 March 2017
     - global gnirs_2017 use
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = gnirs_2017a # Mod on 12/03/2017
    #targets0 = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        run(path0+target+'/', out_script=path0+'run_cleanir.'+target+'.sh')

#enddef
