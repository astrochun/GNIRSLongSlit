"""
QA_plot
=======

Plot a set of FITS images for easier visualization of problems
"""

import sys, os

from chun_codes import systime
from chun_codes import match_nosort_str

from os.path import exists
from astropy.io import ascii as asc
from astropy.io import fits

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import glob

from astropy.table import Table
from astropy import log
from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()
import aplpy

bbox_props = dict(boxstyle="square,pad=0.15", fc="w", alpha=0.5, ec="none")

def main(file_list, path0='', out_pdf='', silent=False, verbose=True):

    '''
    main() function to read in each FITS image and display it on a zscale
    using aplpy

    Parameters
    ----------
    file_list : str
      Filename of ASCII file containing images to be plotted. Use 'all.lis'
      from create_list() since this will include all that will be processed
      Do NOT include full path

    path0 : str
      Full path to where output PDF and FITS file are located. Must end
      with a '/'

    out_pdf : str
      Filename for output PDF. Do NOT include full path

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
	  
    Returns
    -------
    multi-page PDF plot

    Notes
    -----
    Created by Chun Ly, 7 March 2017
    '''

    if silent == False: log.info('### Begin main : '+systime())

    if silent == False: log.info('## Reading : '+path0+file_list)    
    files   = np.loadtxt(path0+file_list, dtype=type(str)).tolist()
    n_files = len(files)

    if out_pdf == '':
        out_pdf = path0+'QA_plot.pdf'
    else:
        out_pdf = path0+out_pdf

    hdr_info_file = path0+'hdr_info.tbl'
    if exists(hdr_info_file):
        if silent == False: log.info('## Reading : '+hdr_info_file)
        tab0 = asc.read(hdr_info_file, format='fixed_width_two_line')
        idx1, idx2 = match_nosort_str(files, tab0['filename'])
        tab0 = tab0[idx2]
    else:
        if silent == False: log.warn('## File not found : '+hdr_info_file)
            
    pp = PdfPages(out_pdf)

    for nn in xrange(n_files):
        if silent == False: log.info('## Reading : '+files[nn])    
        # h0  = fits.getheader(path0+files[nn], 0)
        # im0 = fits.getdata(path0+files[nn], 1)
        hdu0 = fits.open(path0+files[nn])
        im0  = hdu0[1].data
        hdu0[1].header = hdu0[0].header # Copy WCS header over

        gc  = aplpy.FITSFigure(hdu0, hdu=1, figsize=(8,8))

        z1, z2 = zscale.get_limits(im0)
        gc.show_grayscale(invert=True, vmin=z1, vmax=z2)
        gc.set_tick_color('black')

        gc.set_tick_yspacing('auto')
        gc.ticks.set_yspacing(1/60.0) # Every 1 arcmin in Dec
        gc.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm')

        txt0 = files[nn]
        if exists(hdr_info_file):
            txt0 += '\n'
            tmp   = tab0[nn]
            txt0 += 'Date Label : '+tmp['datelabel']+'\n'
            txt0 += 'UT Date : '+tmp['UT_date']+'\n'
            txt0 += tmp['object']+'\n'
            txt0 += 'EXPTIME=%.1f AIRMASS=%.3f \n' % (tmp['exptime'],tmp['airmass'])
            txt0 += '%s %.3f %s' % (tmp['grating'],tmp['gratwave'],tmp['filter2'])

        gc.add_label(0.975, 0.115, txt0, color='red', relative=True, bbox=bbox_props,
                     ha='right', va='bottom', weight='medium', size='medium')

        gc.savefig(pp, format='pdf')

    if silent == False: log.info('## Reading : '+out_pdf)
    pp.close()

    if silent == False: log.info('### End main : '+systime())
#enddef

