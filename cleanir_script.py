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

import glog # + on 22/01/2018

# + on 12/03/2017
from . import gnirs_2017a #targets0
import dir_check # + on 15/05/2017

co     = __file__
co_dir = os.path.dirname(co)+'/'
cmd0   = co_dir+'cleanir.py'
        
def run(path0, clean_file='', out_script='', silent=False, verbose=True,
        overwrite=False):

    '''
    Create a .sh script to run cleanir.py for a set of files

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    overwrite : boolean
      Overwrite files if they exists. Default: False

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 7 March 2017
    Modified by Chun Ly, 15 May 2017
     - Call dir_check.main() to handle multiple date directories
    Modified by Chun Ly, 30 May 2017
     - Added overwrite option. Default is to not overwrite .sh files
     - Fix bug when only one file is found
    Modified by Chun Ly, 18 June 2017
     - Fix to work with no date directory
    Modified by Chun Ly, 22 January 2018
     - Import glog and call for stdout and ASCII logging
     - Pass mylogger to dir_check.main()
    '''

    logfile  = path0+'cleanir_script.log'
    mylogger = glog.log0(logfile)._get_logger()
    
    if silent == False: mylogger.info('### Begin run : '+systime())

    if clean_file == '': clean_file = 'clean.lis'

    # + on 15/05/2017
    dir_list, list_path = dir_check.main(path0, mylogger=mylogger, silent=silent,
                                         verbose=verbose)

    # Mod on 15/05/2017
    for date,path in zip(dir_list,list_path):
        clean_file0 = path+clean_file
        if not exists(clean_file0):
            mylogger.warn('File does not exist!!!')
            mylogger.warn(clean_file0)
        else:
            if silent == False: mylogger.info('Reading : '+clean_file0)
            files = np.loadtxt(clean_file0, dtype=type(str)).tolist()
            if type(files) == str: files = [files] # Bug fix. Mod on 30/05/2017

            out_script0 = path+'run_cleanir.sh' if out_script == '' \
                          else out_script

            if date != '': # Mod on 18/06/2017
                out_script0 = out_script0.replace('.sh', '.'+date+'.sh')

            # Mod on 30/05/2017
            if overwrite == False and exists(out_script0):
                log.warn('## File found!!! : '+out_script0)
                log.warn('## Will not overwrite!!!')
            else:
                if silent == False:
                    stat0 = 'Overwriting' if exists(out_script0) else 'Writing'
                    mylogger.info(stat0+' : '+out_script0)
                f = open(out_script0, 'w')
                for ii in xrange(len(files)):
                    cmd1 = cmd0+' -afqo '+path+'c'+files[ii]+' '+path+files[ii]
                    f.write(cmd1+'\n')
                f.close()
            #endelse
        #endelse
    #endfor

    if silent == False: mylogger.info('End run : '+systime())
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
