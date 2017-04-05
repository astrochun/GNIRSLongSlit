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

# + on 04/04/2017 to handle bright stars
from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()

from ccdproc import cosmicray_lacosmic # + on 01/04/2017

from . import gnirs_2017a
import dir_check

# + on 30/03/2017
co_filename = __file__
bpm_file    = os.path.dirname(co_filename)+'/gnirsn_2012dec05_bpm.fits.gz'
bpm_data    = fits.getdata(bpm_file)
bpmy, bpmx  = np.where(bpm_data == 1)

size2d = u.Quantity((560, 770), u.pixel) # u.Quantity((150, 770), u.pixel)
pos0   = (503, 437)

def get_slit_trace(infile): #, xmin, xmax):
    # Mod on 04/04/2017, aesthetics, fix bug
    # Mod on 04/04/2017, handle CRs affecting trace

    im0, hdr0 = fits.getdata(infile, header=True)

    # + on 04/04/2017
    im0_clean = cosmicray_lacosmic(im0, sigclip=10)
    im0_clean = im0_clean[0]

    y_med0 = np.median(im0_clean, axis=1) # Mod on 04/04/2017
    cen0   = (np.where(y_med0 == np.max(y_med0))[0])[0]

    dy = 20 # + on 04/04/2017
    im0_crop = im0_clean[cen0-dy:cen0+dy,:] # Mod on 04/04/2017

    y_idx, x_idx = np.where(im0_crop >= 0.25*np.max(im0_crop))
    xmin, xmax   = np.min(x_idx), np.max(x_idx)

    dx = 2
    x0 = np.arange(xmin, xmax, dx)

    y0_lo = np.zeros(len(x0))
    y0_hi = np.zeros(len(x0))

    for xx in xrange(len(x0)):
        im0_crop = im0_clean[cen0-dy:cen0+dy,x0[xx]:x0[xx]+dx] # Mod on 04/04/2017
        y_med    = np.median(im0_crop, axis=1)
        edge_idx = np.where(y_med >= 0.1*np.max(y_med))[0]
        if len(edge_idx) > 2:
            y0_lo[xx] = cen0-dy+edge_idx[0]
            y0_hi[xx] = cen0-dy+edge_idx[-1]
        else:
            y0_lo[xx] = np.nan
            y0_hi[xx] = np.nan

    return x0+dx/2.0, y0_lo, y0_hi
#enddef

def mask_bad_pixels(im0):
    # + on 30/03/2017

    im0[bpmy,bpmx] = np.nan
    return im0

def find_gnirs_window_mean(infile):
    # + on 04/04/2017
    # Using a mean approach to get the GNIRS window

    log.info('## Reading : '+infile)
    im0  = fits.getdata(infile)
    hdr0 = fits.getheader(infile, ext=0)

    im0_clean = cosmicray_lacosmic(im0, sigclip=10)
    im0_clean = im0_clean[0]

    im0_mask = mask_bad_pixels(im0_clean)

    mean_y = np.nanmean(im0_mask, axis=1)
    mean_x = np.nanmean(im0_mask, axis=0)

    i_y = np.where(mean_y > 0)[0]
    i_x = np.where(mean_x > 0)[0]

    # + on 01/04/2017
    y_min, y_max = np.min(i_y), np.max(i_y)
    x_min, x_max = np.min(i_x), np.max(i_x)
    x_cen, y_cen = (x_max+x_min)/2.0, (y_max+y_min)/2.0 # Later mod on 04/04/2017
    #x_cen, y_cen = np.average(i_x), np.average(i_y)

    # + on 01/04/2017
    info0  = 'x_min=%i, x_max=%i, y_min=%i, y_max=%i ' % \
             (x_min, x_max, y_min, y_max)
    info0 += 'x_cen=%.2f, y_cen=%.2f' % (x_cen, y_cen)
    log.info(info0)
    return x_min, x_max, y_min, y_max, x_cen, y_cen
#enddef

def find_gnirs_window(infile):
    # + on 30/03/2017
    # Mod on 01/04/2017 to figure out way to find GNIRS window

    # Mod on 01/04/2017
    log.info('## Reading : '+infile)
    im0  = fits.getdata(infile)
    hdr0 = fits.getheader(infile, ext=0)

    # Mod on 01/04/2017
    im0_clean = cosmicray_lacosmic(im0, sigclip=10)
    im0_clean = im0_clean[0]

    im0_mask = mask_bad_pixels(im0_clean)

    # + on 01/04/2017
    if hdr0['NDAVGS'] == 1:  v_min = 10
    if hdr0['NDAVGS'] == 16: v_min = 100

    i_y, i_x = np.where((im0_mask > v_min) & np.isfinite(im0_mask))
    print len(i_y)
    med0 = np.median(im0_mask[i_y,i_x])

    # + on 01/04/2017
    y_min, y_max = np.min(i_y), np.max(i_y)
    x_min, x_max = np.min(i_x), np.max(i_x)
    x_cen, y_cen = np.average(i_x), np.average(i_y)

    # + on 01/04/2017
    info0  = '## med0=%.2f, x_min=%i, x_max=%i, y_min=%i, y_max=%i ' % \
             (med0, x_min, x_max, y_min, y_max)
    info0 += 'x_cen=%.2f, y_cen=%.2f' % (x_cen, y_cen)
    log.info(info0)
    return med0, x_min, x_max, y_min, y_max, x_cen, y_cen
#enddef

def find_star(infile, pos=pos0, find_size2d=size2d):
    ## Mod on 24/02/2017 to handle CRs (quick fix)
    ## Mod on 01/04/2017 to handle CRs and bad pixels using cosmicrays_lacosmic
    ## Mod on 04/04/2017 to pass pos and find_size2d keywords

    im0, hdr0 = fits.getdata(infile, header=True)

    # + on 01/04/2017
    im0_clean = cosmicray_lacosmic(im0, sigclip=10)

    # Mod on 04/04/2017
    # find_size2d = size2d #u.Quantity((25, 770), u.pixel)
    cutout = Cutout2D(im0_clean[0], pos, find_size2d, mode='partial',
                      fill_value=np.nan)
    cutout = cutout.data

    peak = np.where(cutout == np.max(cutout))
    ycen, xcen = peak[0][0], peak[1][0]
    print xcen, ycen

    xcen += pos[0]-find_size2d[1].value/2.0
    ycen += pos[1]-find_size2d[0].value/2.0
    print xcen, ycen
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
    Modified by Chun Ly, 01 April 2017
     - Handle CRs and bad pixels using cosmicrays_lacosmic
    Modified by Chun Ly, 04 April 2017
     - Use find_gnirs_window_mean to find center
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

        # + on 04/04/2017
        win_ref_idx  = [tt for tt in xrange(len(tab0)) if
                        (tab0['QA'][tt] == 'N/A') and ('Acq' in tab0['slit'][tt]) and
                        ('HIP' not in tab0['object'][tt]) and
                        ('HD' not in tab0['object'][tt])]
        if len(win_ref_idx) > 0:
            win_ref_file = path + tab0['filename'][win_ref_idx[0]]
            log.info('## Reference image for finding GNIRS window : '+win_ref_file)

            x_min, x_max, y_min, y_max, \
                x_cen, y_cen = find_gnirs_window_mean(win_ref_file)
            pos_cen  = (x_cen, y_cen)
            new_size = u.Quantity((y_max-y_min, x_max-x_min), u.pixel)
            print 'pos_cen : ', pos_cen
            print 'new_size : ', new_size

        for ii in xrange(len(ID0)):
            t_idx = [tt for tt in xrange(len(tab0)) if
                     (tab0['object'][tt] == ID0[ii] and
                      tab0['QA'][tt] == 'N/A')]

            t_files = [path+a for a in tab0['filename'][t_idx]]
            ncols = 2.0
            nrows = np.ceil(len(t_idx)/ncols)
            ncols, nrows = np.int(ncols), np.int(nrows)

            fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols)

            #med0, x_min, x_max, y_min, \
            #    y_max, x_cen, y_cen = find_gnirs_window(t_files[1])

            # Later + on 24/03/2017 | Mod on 04/04/2017
            xcen, ycen = find_star(t_files[-1], pos=pos_cen, find_size2d=new_size)
            # Fix to get relative coordinate for Cutout2D image
            #xcen -= pos_cen[0]-new_size[1].value/2.0
            #ycen -= pos_cen[1]-new_size[0].value/2.0

            slit_x0, slit_y0_lo, slit_y0_hi = get_slit_trace(t_files[0]) #, x_min, x_max)
            # Adjust values for offset that is applied
            # Bug: Mod on 04/04/2017 to get proper coordinate
            slit_x0    -= np.int64(pos_cen[0]-size2d[1].value/2.0)
            slit_y0_lo -= pos_cen[1]-size2d[0].value/2.0
            slit_y0_hi -= pos_cen[1]-size2d[0].value/2.0

            for jj in xrange(len(t_idx)):
                jj_idx = t_idx[jj]

                im0, hdr0 = fits.getdata(t_files[jj], header=True)

                # + 01/04/2017
                im0_clean = cosmicray_lacosmic(im0, sigclip=10)

                cutout = Cutout2D(im0_clean[0], pos_cen, size2d,
                                  mode='partial', fill_value=np.nan)

                t_col, t_row = jj % ncols, jj / ncols

                # Mod on 04/04/2017 to handle bright and faint stars
                max0 = np.max(cutout.data)
                if max0 > 50000:
                    z1, z2 = zscale.get_limits(cutout.data)
                else:
                    z1, z2 = 0.0, 0.5*max0

                norm = ImageNormalize(vmin=z1, vmax=z2)
                t_ax = ax_arr[t_row,t_col]
                t_ax.imshow(cutout.data, cmap='Greys', origin='lower',
                            norm=norm)
                #aplpy.FITSFigure(cutout)

                # Draw trace of slit
                t_ax.plot(slit_x0, slit_y0_lo, 'r-')
                t_ax.plot(slit_x0, slit_y0_hi, 'r-')

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
                norm2 = ImageNormalize(vmin=z1, vmax=z2)
                axins.imshow(cutout.data, cmap='Greys', origin='lower',
                             norm=norm2)

                # Draw trace of slit
                axins.plot(slit_x0, slit_y0_lo, 'r-')
                axins.plot(slit_x0, slit_y0_hi, 'r-')

                # Mod on 04/04/2017 to get Cutout2d coordinates
                c_xcen = xcen - (pos_cen[0]-size2d[1].value/2.0)
                c_ycen = ycen - (pos_cen[1]-size2d[0].value/2.0)
                x1, x2, y1, y2 = c_xcen-20, c_xcen+20, c_ycen-20, c_ycen+20
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
