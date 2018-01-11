"""
wave_cal_script
===============

Generate a PyRAF script to run interactive wavelength calibration
"""

import sys, os
from chun_codes import systime
from os.path import exists
import numpy as np
from astropy import log

from check_path import main as check_path

n_sp_pix = 1022
yes, no  = 'yes', 'no'

from astropy.io import fits

import glog # + on 09/01/2018

co_dirname = os.path.dirname(__file__)
OH_file = co_dirname+'/rousselot2000.dat'

def main(rawdir, line_source='', mylogger=None, silent=False, verbose=True):

    '''
    Main function for wave_cal_script

    Parameters
    ----------
    rawdir : str
      Path to raw files. Must end in a '/'

    line_source : str
      Type of lines to use for wavelength calibration. Option is
      either 'arc' or 'OH

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 7-8 October 2017
    Modified by Chun Ly, 8 November 2017
     - Specify lampspec and outspec
    Modified by Chun Ly, 16 November 2017
     - Handle line_source == 'OH'
     - Generalized arrays: arc -> frame
     - Specify GNIRS log file
    Modified by Chun Ly, 20 November 2017
     - Bug fix: prefix needed for [frames] for line_source == OH
    Modified by Chun Ly,  9 January 2018
     - Import glog and call for stdout and ASCII logging
    Modified by Chun Ly,  9 January 2018
     - Allow mylogger keyword
    '''

    # Mod on 10/01/2018
    if type(mylogger) == type(None):
        logfile  = rawdir+'QA_wave_cal.log'
        mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin main : '+systime())

    timestamp = systime().replace(':','.')
    logfile   = rawdir+'gnirs_'+timestamp+'.log'

    rawdir = check_path(rawdir)

    if line_source == '':
        mylogger.warn("Must specify line_source keyword: ")
        mylogger.warn("line_source='arc' or line_source='OH'")
        mylogger.warn("Exiting!!!")
        return
    else:
        out_script = rawdir + 'wave_cal_'+line_source+'.py'

    if silent == False: mylogger.info('Writing : '+out_script)

    f0 = open(out_script, 'w')

    line0 = ['iraf.gemini(_doprint=0)', 'iraf.gemini.gnirs(_doprint=0)',
             'iraf.gemini.unlearn()', 'iraf.gemini.gemtools.unlearn()',
             'iraf.gemini.gnirs.unlearn()', 'iraf.set(stdimage="imt4096")',
             'iraf.gemini.nsheaders("gnirs")']
    #'iraf.gemini.gnirs.logfile = "%s"' % logfile] # + on 16/11/2017

    f0.writelines("\n".join(line0)+"\n")

    if line_source == 'arc': frame_list = rawdir+'arc.lis'
    if line_source == 'OH':  frame_list = rawdir+'obj.OH.lis'

    frames    = np.loadtxt(frame_list, dtype=type(str))
    if line_source == 'OH':
        frames = ['rbnc'+file0 for file0 in frames]

    frame_hdr = fits.getheader(rawdir+frames[0])

    crpix = n_sp_pix / 2.0
    crval = frame_hdr['gratwave'] * 1e4 # in Angstroms
    if frame_hdr['FILTER2'] == 'X_G0518':
        cdelt = -0.094*1e4/n_sp_pix
    if frame_hdr['FILTER2'] == 'J_G0517':
        cdelt = -0.113*1e4/n_sp_pix
    mylogger.info('## CRVAL : %.1f ' % crval)
    mylogger.info('## CDELT : %.1f  CRPIX : %.1f' % (cdelt,crpix))

    line1 = ['crval = %f' % crval, 'crpix = %f' % crpix, 'cdelt = %f' % cdelt]
    f0.writelines("\n".join(line1)+"\n")

    if line_source == 'arc':
        coordlist = 'gnirs$data/argon.dat'
        database  = 'database/'
    if line_source == 'OH':
        coordlist = OH_file
        database  = 'database_OH/'

    # Mod on 08/11/2017
    line2 = ["coordlist = '%s'" % coordlist,
             "database  = '%s'" % database,
             "lampspec  = '%s_stack.fits'"  % line_source,
             "outspec   = 'w%s_stack.fits'" % line_source,
             "logfile   = '%s'" % logfile] # + on 16/11/2017

    f0.writelines("\n".join(line2)+"\n")

    # Mod on 08/11/2017, 16/11/2017
    cmd = "iraf.gnirs.nswavelength(lampspec, outprefix='',"+\
          "outspectra=outspec, crval=crval, cdelt=cdelt, crpix=crpix, "+\
          "coordlist=coordlist, database=database, fl_inter='yes',"+\
          "cradius=20, threshold=50.0, order=2, logfile=logfile)"

    f0.write(cmd+'\n')
    f0.close()

    if silent == False: mylogger.info('### End main : '+systime())
#enddef
