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

def main(path0, silent=False, verbose=True):
    '''
    Main function to get list of directories formatted as
    YYYYMMDD

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------
    dir_list : list
     List of date directories in [path0]

    Notes
    -----
    Created by Chun Ly, 23 March 2017
    '''

    if silent == False: log.info('### Begin main : '+systime())

    list0 = [os.path.join(path0,f) for f in os.listdir(path0)]
    dirs = filter(os.path.isdir, list0)

    dir_list = []
    for dd in xrange(len(dirs)):
        t_date = dirs[dd].split('/')[-1]
        chk = is_date(t_date)
        if chk:
            dir_list.append(t_date)

    if silent == False: log.info('### End main : '+systime())
    return dir_list
#enddef

