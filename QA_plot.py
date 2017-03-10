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

from astropy.stats import sigma_clipped_stats # + on 10/03/2017

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
    Created by Chun Ly, 6 March 2017
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
            txt0 += 'EXPTIME=%.1f ' % tmp['exptime']
            txt0 += 'AIRMASS=%.3f \n' % tmp['airmass']
            txt0 += '%s %.3f %s' % (tmp['grating'],tmp['gratwave'],
                                    tmp['filter2'])

        gc.add_label(0.975, 0.115, txt0, color='red', relative=True, ha='right',
                     va='bottom', weight='medium', size='medium',
                     bbox=bbox_props)

        gc.savefig(pp, format='pdf')

    if silent == False: log.info('## Reading : '+out_pdf)
    pp.close()

    if silent == False: log.info('### End main : '+systime())
#enddef

def quadrant_bias_values(hdu, gc0):
    '''
    Compute median, mean, and average and annotate subplots

    Parameters
    ----------
    hdu : an `HDUList` object
      `HDUList` containing all of the header data units in the file.

    gc0 : aplpy FITSFigure
      A pointer to the aplpy.FITSFigure class

    Returns
    -------
    Annotates subplots through gc0

    Notes
    -----
    Created by Chun Ly, 10 March 2017
     - Later added to compute values using a sigma-clipping approach
    '''
    bad = -9.e6

    next = len(hdu)
    if ( next == 1 ):
        sci = 0
    elif ( next < 5 ):
        sci = 1
    else:
        sci = 2

    naxis1, naxis2 = hdu[sci].header['naxis1'], hdu[sci].header['naxis2']

    qxsize = naxis1 / 2 # quadrant x size
    qysize = naxis2 / 2 # quadrant y size

    q1x1,q1x2, q1y1,q1y2 = 0,qxsize,      qysize,naxis2  # quad 1
    q2x1,q2x2, q2y1,q2y2 = qxsize,naxis1, qysize,naxis2  # quad 2
    q3x1,q3x2, q3y1,q3y2 = 0,qxsize,           0,qysize  # quad 3
    q4x1,q4x2, q4y1,q4y2 = qxsize,naxis1,      0,qysize  # quad 4

    im0 = np.array(hdu[sci].data)

    quads = im0.copy()
    quad1 = quads[qysize:naxis2,      0:qxsize]
    quad2 = quads[qysize:naxis2, qxsize:naxis1]
    quad3 = quads[    0:qysize,       0:qxsize]
    quad4 = quads[    0:qysize,  qxsize:naxis1]

    str0 = 'CLEANIR Methodology\n'
    str1 = 'Sigma-clipped Methodology\n'

    for qq,quad0 in zip(range(1,5),[quad1,quad2,quad3,quad4]):
        mean0 = np.mean(quad0)
        med0  = np.median(quad0)
        sig0  = np.std(quad0)
        str0 += 'Q%i mean=%.3f median=%.3f sigma=%.3f' % (qq,mean0,med0,sig0)
        if qq != 4: str0 += '\n'

        # Later + on 10/03/2017
        mean0, med0, sig0 = sigma_clipped_stats(quad0, sigma=2.0, iters=10)
        str1 += 'Q%i mean=%.3f median=%.3f sigma=%.3f' % (qq,mean0,med0,sig0)
        if qq != 4: str1 += '\n'

    gc0.add_label(0.99,0.060, str0, color='red', relative=True, ha='right',
                  va='bottom', weight='medium', size='medium', bbox=bbox_props)
    # Later + on 10/03/2017
    gc0.add_label(0.99,0.200, str1, color='red', relative=True, ha='right',
                  va='bottom', weight='medium', size='medium', bbox=bbox_props)

#enddef

def clean_QA(path0='', out_pdf='', silent=False, verbose=True):
    '''
    Visually compare raw and cleanir-fixed FITS images

    Parameters
    ----------
    path0 : str
      Full path to where output PDF and inputs FITS file are located. Must end
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
    Created by Chun Ly, 8 March 2017
    Modified by Chun Ly, 10 March 2017
     - Call quadrant_bias_values()
    '''

    if silent == False: log.info('### Begin clean_QA : '+systime())

    if out_pdf == '':
        out_pdf = path0+'QA_plot.pdf'
    else:
        out_pdf = path0+out_pdf

    files   = glob.glob(path0+'cN*fits')
    n_files = len(files)

    pp = PdfPages(out_pdf)

    for nn in xrange(n_files):
        if silent == False: log.info('## Reading : '+files[nn])
        orig_file = files[nn].replace('cN','N')

        im1 = fits.getdata(orig_file)
        im2 = fits.getdata(files[nn])
        #hdu1 = fits.open(orig_file)
        #im1  = hdu1[1].data
        #hdu1[1].header = hdu1[0].header # Copy WCS header over

        #hdu2 = fits.open(files[nn])
        #im2  = hdu2[1].data
        #hdu2[1].header = hdu2[0].header # Copy WCS header over

        fig = plt.figure(figsize=(16,8))

        gc1 = aplpy.FITSFigure(orig_file, figure=fig,
                               subplot=[0.05,0.055,0.47,0.94])
        z1, z2 = zscale.get_limits(im1)
        gc1.show_grayscale(invert=True, vmin=z1, vmax=z2)
        gc1.set_tick_color('black')
        gc1.set_tick_yspacing('auto')
        #gc1.ticks.set_yspacing(1/60.0) # Every 1 arcmin in Dec
        #gc1.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm')
        gc1.add_label(0.025, 0.975, orig_file, color='red', relative=True,
                      ha='left', va='top', weight='medium', size='medium',
                      bbox=bbox_props)

        quadrant_bias_values(fits.open(orig_file), gc1) # + on 10/03/2017

        gc2 = aplpy.FITSFigure(files[nn], figure=fig,
                               subplot=[0.525,0.055,0.47,0.94])
        z1, z2 = zscale.get_limits(im2)
        gc2.show_grayscale(invert=True, vmin=z1, vmax=z2)
        gc2.set_tick_color('black')
        gc2.set_tick_yspacing('auto')
        gc2.hide_ytick_labels()
        gc2.hide_yaxis_label()
        gc2.add_label(0.025, 0.975, files[nn], color='red', relative=True,
                      ha='left', va='top', weight='medium', size='medium',
                      bbox=bbox_props)

        quadrant_bias_values(fits.open(files[nn]), gc2) # + on 10/03/2017

        #gc2.ticks.set_yspacing(1/60.0) # Every 1 arcmin in Dec
        #gc2.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm')

        fig.savefig(pp, format='pdf')

    if silent == False: log.info('## Reading : '+out_pdf)
    pp.close()

    if silent == False: log.info('### End clean_QA : '+systime())
#enddef

def zcalbase_gal_gemini_2017a_raw():
    '''
    Function to run main() on each set of GNIRS 2017A observation set
    to obtain PDF plots of each *raw* FITS image for visual examination

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 6 March 2017
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        main('all.lis', path0=path0+target+'/', out_pdf='QA_plot.raw.pdf')

#enddef

def zcalbase_gal_gemini_2017a_cleanir():
    '''
    Function to run main() on each set of GNIRS 2017A observation set
    to obtain PDF plots of each *raw* FITS image for visual examination

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 6 March 2017
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

    for target in targets0:
        clean_QA(path0=path0+target+'/', out_pdf='QA_plot.clean.pdf')

#enddef
