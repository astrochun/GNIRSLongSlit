"""
check_path
==========

Fix paths to avoid problem with '/' requirement
"""

import sys, os

from chun_codes import systime

from os.path import exists

def main(path0, silent=True, verbose=False):

    '''
    Main function for checking paths and adding '/'

    Parameters
    ----------
    path0 : str
      Full path.

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 20 September 2017
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if path0[-1] != '/': path0 = path0 + '/'

    if silent == False: log.info('### End main : '+systime())
    return path0
#enddef

