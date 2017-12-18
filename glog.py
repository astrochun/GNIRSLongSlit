"""
glog
====

Routine to log information to stdout and ASCII file

To execute:
from GNIRSLongSlit import glog

logfile = '/path/mylog.log'
mylogger = glog.log0(logfile)._get_logger()
"""

import logging, sys

# LOG_FILENAME = '/Users/cly/test.log'

formatter = logging.Formatter('%(asctime)s - %(module)12s.%(funcName)20s - %(levelname)s: %(message)s')

# set up logging to STDOUT for all levels INFO and higher
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
sh.setFormatter(formatter)
#sh.handler_set = True

class log0:
    '''
    Main class to log information to stdout and ASCII file

    To execute:
    logfile = '/path/mylog.log'
    mylogger = log0(logfile)._get_logger()

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

    Modified by Chun Ly, 8 December 2017
     - Switch from function to class to avoid repeated stdout entries due to
       calling of handlers
    Modified by Chun Ly, 14 December 2017
     - Fix class instance call to write to separate files with multiple calls
    Modified by Chun Ly, 17 December 2017
     - Documentation update
    '''

    def __init__(self,file):
        self.LOG_FILENAME = file
        self._log = self._get_logger()

    def _get_logger(self):
        loglevel = logging.INFO
        log = logging.getLogger(self.LOG_FILENAME) # + Mod on 14/12/2017
        if not getattr(log, 'handler_set', None):
            log.setLevel(logging.INFO)
            sh = logging.StreamHandler()
            sh.setFormatter(formatter)
            log.addHandler(sh)

            fh = logging.FileHandler(self.LOG_FILENAME)
            fh.setLevel(logging.INFO)
            fh.setFormatter(formatter)
            log.addHandler(fh)

            log.setLevel(loglevel)
            log.handler_set = True
        return log

    #debug = mylogger.debug
    #info = mylogger.info
    #warning = mylogger.warning
    #error = mylogger.error
    #critical = mylogger.critical

    #return mylogger
#enddef
