"""
glog
====

Routine to log information to stdout and ASCII file
"""

import logging, sys

# LOG_FILENAME = '/Users/cly/test.log'

formatter = logging.Formatter('%(asctime)s - %(module)12s.%(funcName)20s - %(levelname)s: %(message)s')

# set up logging to STDOUT for all levels INFO and higher
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)
#sh.handler_set = True


def log0(LOG_FILENAME):
    '''
    Main function to log information to stdout and ASCII file

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
    Created by Chun Ly, 7 December 2017
    '''

    fh = logging.FileHandler(LOG_FILENAME)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)

    # create Logging object
    mylogger = logging.getLogger('MyLogger')
    mylogger.setLevel(logging.DEBUG)
    mylogger.addHandler(sh)    # enabled: stdout
    mylogger.addHandler(fh)    # enabled: file
    mylogger.propagate = False

    #debug = mylogger.debug
    #info = mylogger.info
    #warning = mylogger.warning
    #error = mylogger.error
    #critical = mylogger.critical

    return mylogger
#enddef
