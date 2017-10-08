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

co_dirname = os.path.dirname(__file__)
OH_file = co_dirname+'/rousselot2000.dat'

def main(rawdir, line_source='', silent=False, verbose=True):

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
    '''

    if silent == False: log.info('### Begin main : '+systime())

    rawdir = check_path(rawdir)

    arc_list  = rawdir+'arc.lis'

    if line_source == '':
        log.warn("## Must specify line_source keyword: ")
        log.warn("## line_source='arc' or line_source='OH'")
        log.warn("## Exiting!!!")
        return
    else:
        out_script = rawdir + 'wave_cal_'+line_source+'.py'

    if silent == False: log.info('### Writing : '+out_script)

    f0 = open(out_script, 'w')

    line0 = ['iraf.gemini(_doprint=0)', 'iraf.gemini.gnirs(_doprint=0)',
             'iraf.gemini.unlearn()', 'iraf.gemini.gemtools.unlearn()',
             'iraf.gemini.gnirs.unlearn()', 'iraf.set(stdimage="imt4096")',
             'iraf.gemini.nsheaders("gnirs")']

    f0.writelines("\n".join(line0)+"\n")

    arcs = np.loadtxt(arc_list, dtype=type(str))
    arc_hdr = fits.getheader(rawdir+arcs[0])

    crpix    = n_sp_pix / 2.0
    crval    = arc_hdr['gratwave'] * 1e4 # in Angstroms
    if arc_hdr['FILTER2'] == 'X_G0518':
        cdelt = -0.094*1e4/n_sp_pix
    if arc_hdr['FILTER2'] == 'J_G0517':
        cdelt = -0.113*1e4/n_sp_pix
    log.info('## CRVAL : %.1f ' % crval)
    log.info('## CDELT : %.1f  CRPIX : %.1f' % (cdelt,crpix))

    line1 = ['crval = %f' % crval, 'crpix = %f' % crpix, 'cdelt = %f' % cdelt]
    f0.writelines("\n".join(line1)+"\n")

    if line_source == 'arc':
        coordlist = 'gnirs$data/argon.dat'
        database='database/'
    if line_source == 'OH':
        coordlist = OH_file
        database='database_OH/'
        
    line2 = ['coordlist = %s' % coordlist, 'database = %s' % database]
    f0.writelines("\n".join(line2)+"\n")

    cmd = "iraf.gnirs.nswavelength('rnc@arc.lis', outprefix='',"+\
          "outspectra='wrnc@arc.lis', crval=crval, cdelt=cdelt, crpix=crpix"+\
          "coordlist=coordlist, database=database, fl_inter='yes',"+\
          "cradius=20, threshold=50.0, order=2)"

    f0.write(cmd)
    f0.close()

    if silent == False: log.info('### End main : '+systime())
#enddef

