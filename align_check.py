"""
align_check
===========

Illustrate alignment using alignment FITS files
"""

import sys, os

from chun_codes import systime
import numpy as np

from os.path import exists
from astropy.io import ascii as asc
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pylab import subplots_adjust
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import glob

import astropy.units as u
from astropy import log
from astropy.nddata import Cutout2D
from astropy.visualization.mpl_normalize import ImageNormalize

from . import gnirs_2017a
import dir_check

size2d = u.Quantity((560, 770), u.pixel) # u.Quantity((150, 770), u.pixel)
pos0   = (503, 437)

def find_star(infile):
    im0, hdr0 = fits.getdata(infile, header=True)

    cutout = Cutout2D(im0, pos0, size2d, mode='partial',
                      fill_value=np.nan)
    cutout = cutout.data

    peak = np.where(cutout == np.max(cutout))
    ycen, xcen = peak[0], peak[1]
    return xcen, ycen
#enddef

def main(path0, out_pdf='', silent=False, verbose=True):

    '''
    Main function to generate PDF illustrating alignment on target

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

    Notes
    -----
    Created by Chun Ly, 24 March 2017
    '''
    
    if silent == False: log.info('### Begin main : '+systime())
    
    dir_list, list_path = dir_check.main(path0, silent=silent, verbose=verbose)

    out_pdf_default = out_pdf
    for path in list_path:
        infile = path + 'hdr_info.QA.tbl'
        if not exists(infile):
            log.warning('### File does not exist : '+infile)
            log.warning('### Exiting!!! '+systime())
            return

        out_pdf = path+'align_check.pdf' if out_pdf == '' else path+out_pdf
        pp = PdfPages(out_pdf)

        if silent == False: log.info('### Reading: '+infile)
        tab0 = asc.read(infile, format='fixed_width_two_line')

        align = [ii for ii in xrange(len(tab0)) if tab0['QA'][ii] == 'N/A']
        if silent == False:
            log.info('## Number of alignment images found : '+str(len(align)))

        ID  = tab0['object'][align]
        ID0 = list(set(ID)) #Unique ID's
        if silent == False:
            log.info('## Sources found : '+', '.join(ID0))

        for ii in xrange(len(ID0)):
            t_idx = [tt for tt in xrange(len(tab0)) if
                     (tab0['object'][tt] == ID0[ii] and
                      tab0['QA'][tt] == 'N/A')]

            t_files = [path+a for a in tab0['filename'][t_idx]]
            ncols = 2.0
            nrows = np.ceil(len(t_idx)/ncols)
            ncols, nrows = np.int(ncols), np.int(nrows)

            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols)

            # Later + on 24/03/2017
            xcen, ycen = find_star(t_files[-1])

            for jj in xrange(len(t_idx)):
                jj_idx = t_idx[jj]

                im0, hdr0 = fits.getdata(t_files[jj], header=True)
                cutout = Cutout2D(im0, pos0, size2d, mode='partial',
                                  fill_value=np.nan)

                t_col, t_row = jj % ncols, jj / ncols

                # med0 = np.absolute(np.median(cutout.data))
                # print med0
                max0 = np.max(cutout.data)
                norm = ImageNormalize(vmin=0.0, vmax=0.1*max0)
                t_ax = ax_arr[t_row,t_col]
                t_ax.imshow(cutout.data, cmap='Greys', origin='lower',
                            norm=norm)
                #aplpy.FITSFigure(cutout)

                t_ax.xaxis.set_ticklabels([])
                t_ax.yaxis.set_ticklabels([])

                fig.suptitle(path, fontsize=14)

                txt0  = tab0['filename'][jj_idx]+'\n'
                txt0 += tab0['datelabel'][jj_idx]+'\n'
                txt0 += tab0['UT_date'][jj_idx]+'\n'
                txt0 += tab0['object'][jj_idx]
                t_ax.annotate(txt0, [0.025, 0.95], xycoords='axes fraction',
                              ha='left', va='top')

                # Plot inset | Later + on 24/03/2017
                axins = zoomed_inset_axes(t_ax, 3, loc=4)
                norm2 = ImageNormalize(vmin=0.0, vmax=0.2*max0)
                axins.imshow(cutout.data, cmap='Greys', origin='lower',
                             norm=norm2)

                x1, x2, y1, y2 = xcen-20, xcen+20, ycen-20, ycen+20
                axins.set_xlim([x1, x2])
                axins.set_ylim([y1, y2])
                axins.xaxis.set_ticklabels([])
                axins.yaxis.set_ticklabels([])
                mark_inset(t_ax, axins, loc1=1, loc2=3, fc="none", ec="b",
                           ls='dotted', lw=0.5)
            #endfor

            subplots_adjust(left=0.02, bottom=0.02, top=0.95, right=0.98,
                            wspace=0.02, hspace=0.02)
            fig.set_size_inches(11,8)
            fig.savefig(pp, format='pdf')
        #endfor

        pp.close()
        out_pdf = out_pdf_default
        
    if silent == False: log.info('### End main : '+systime())
#enddef

def zcalbase_gal_gemini_2017a():
    '''
    Function to run main() on each GNIRS 2017A observation set
    to generate visualization of slit alignment on bright offset star

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Created by Chun Ly, 24 March 2017
    '''

    path0 = '/Users/cly/data/Observing/Gemini/Data/'

    targets0 = gnirs_2017a

    for target in targets0:
        main(path0+target+'/')
#enddef
