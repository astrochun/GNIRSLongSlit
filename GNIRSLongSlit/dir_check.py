"""
dir_check
=========

Provide description for code here.
"""

import sys, os

from chun_codes import systime
from os.path import exists
import glob

from astropy import log
from dateutil.parser import parse

def is_date(string):
    '''
    Check if a string is a date
    '''

    try: 
        parse(string)
        return True
    except ValueError:
        return False

def main(path0, mylogger=None, silent=False, verbose=True):
    '''
    Main function to get list of directories formatted as
    YYYYMMDD

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
    dir_list : list
      List of date directories in [path0]

    list_path : list
      List containing full path for date directories

    Notes
    -----
    Created by Chun Ly, 23 March 2017
    Modified by Chun Ly, 18 June 2017
     - Change dir_list output if no directories is available
    Modified by Chun Ly, 10 December 2017
    - Implement glog logging, allow mylogger keyword input
    '''

    # + on 10/12/2017
    if type(mylogger) == type(None):
        mylog, clog = 0, log
    else:
        mylog, clog = 1, mylogger

    if silent == False: clog.info('### Begin main : '+systime())

    list0 = [os.path.join(path0,f) for f in os.listdir(path0)]
    dirs = filter(os.path.isdir, list0)

    dir_list = []
    for dd in xrange(len(dirs)):
        t_date = dirs[dd].split('/')[-1]
        chk = is_date(t_date)
        if chk:
            dir_list.append(t_date)

    # Later + on 23/03/2017
    if len(dir_list) == 0:
        if silent == False:
            clog.info('No dir found') # Mod on 10/12/2017

        dir_list = [''] # Mod on 18/07/2017
        list_path = [path0]
    else:
        if silent == False:
            # Mod on 10/12/2017
            txt0 = 'The following date dir found: '+', '.join(dir_list)
            if not mylog: txt0 = '## '+txt0
            clog.info(txt0)

        list_path = [path0+a+'/' for a in dir_list]

    if silent == False: clog.info('### End main : '+systime()) # Mod on 10/12/2017
    return dir_list, list_path
#enddef

