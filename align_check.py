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

# + on 05/04/2017 to group index together
from itertools import groupby
from operator import itemgetter

# For 2-D gaussian fitting | + on 05/04/2017
# Mod on 26/04/2017 to use chun_codes instead
from chun_codes import gauss2d
#MMTtools.mmtcam import gauss2d
import scipy.optimize as opt
f_s = 2*np.sqrt(2*np.log(2)) # sigma-FWHM conversion

# + on 06/04/2017
from scipy.interpolate import interp1d

from . import gnirs_2017a
import dir_check

# + on 30/03/2017
co_filename = __file__
bpm_file    = os.path.dirname(co_filename)+'/gnirsn_2012dec05_bpm.fits.gz'
bpm_data    = fits.getdata(bpm_file)
bpmy, bpmx  = np.where(bpm_data == 1)

size2d = u.Quantity((560, 770), u.pixel) # u.Quantity((150, 770), u.pixel)
pos0   = (503, 437)

def group_index(index0, find_max=False):

    '''
    Function to group sets of index

    Parameters
    ----------
    index0 : list or numpy array
      List or numpy array containing the index

    find_max : bool
      Keyword flag for whether to return the index array with the maximum
      size or the full list

    Returns
    -------
    list0 : list
      List of numpy arrays with each entry a different group index

    list0[m0] : numpy array
      Array containing the largest number of indexes

    Notes
    -----
    Created by Chun Ly, 5 Apri 2017
    '''

    count0 = list()
    list0  = list()
    for k, g in groupby(enumerate(index0), lambda (i, x): i-x):
        temp = np.array(map(itemgetter(1), g))
        list0.append(temp)
        count0.append(len(temp))

    if find_max == False:
        return list0
    else:
        m0 = np.where(count0 == np.max(count0))[0][0]
        return list0[m0]

#enddef

def gauss2d_fit(im0, hdr0, t_ax, c_slit_x0, c_slit_y0_lo, c_slit_y0_hi):
    '''
    Execute 2-D Gaussian fit on image

    Parameters
    ----------
    im0 : numpy array
      Cropped out image for fitting

    hdr0 : FITS header
      FITS header containing WCS information for pixel scale

    t_ax : matplotlib axis pointer
      To indicate which panel to annotate text

    c_slit_x0: array
      Contains arrays of x values

    c_slit_y0_lo: array
      Contains arrays of y values for the base of the slit

    c_slit_y0_hi: array
      Contains arrays of y values for the top of the slit

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 5 April 2017
    Modified by Chun Ly, 6 April 2017
     - Determine how far source is from slit center
    '''
    sigG = 3.0
    max0 = np.nanmax(im0)

    shape0 = im0.shape
    x_sz0, y_sz0 = shape0[0]/2.0, shape0[1]/2.0

    if max0 < 50000:
        pscale = np.sqrt(hdr0['CD1_1']**2 + hdr0['CD2_1']**2)*3600.0*u.arcsec
        pscale = pscale.to(u.arcsec).value

        gx     = np.linspace(0,shape0[0]-1,shape0[0]) - x_sz0
        gy     = np.linspace(0,shape0[1]-1,shape0[1]) - y_sz0
        gx, gy = np.meshgrid(gx, gy)

        im0_re = im0.reshape(shape0[0]*shape0[1])

        m_idx = np.where(im0 == max0)
        x_cen, y_cen = m_idx[1][0]-x_sz0, m_idx[0][0]-y_sz0
        print '## cen : ', x_cen, y_cen
        ini_guess = (max0, x_cen, y_cen, sigG, sigG, 0.0, 0.0)
        bounds = ((     0, -x_sz0, -y_sz0,    0,    0,       0, 0),
                  (np.inf,  x_sz0,  y_sz0, 10.0, 10.0, 2*np.pi, np.inf))
        popt, pcov = opt.curve_fit(gauss2d, (gx, gy), im0_re, ini_guess,
                                   bounds=bounds)

        FWHMx = popt[3] * f_s * pscale
        FWHMy = popt[4] * f_s * pscale

        # Determine how off from slit center | + on 06/04/2017
        f_lo = interp1d(c_slit_x0, c_slit_y0_lo)
        f_hi = interp1d(c_slit_x0, c_slit_y0_hi)
        s_lo = f_lo(popt[1])-y_sz0
        s_hi = f_hi(popt[1])-y_sz0
        slit_off = (popt[2] - (s_lo+s_hi)/2.0) * pscale
        print '## slit_off(arcsec) : ', slit_off
        str0  = 'Slit offset : %.2f"\n' % slit_off

        str0 += 'Centroid: x=%.2f"  y=%.2f"\n' % (popt[1]*pscale, popt[2]*pscale)
        str0 += 'FWHM: x=%.2f"  y=%.2f"\n'   % (FWHMx, FWHMy)
        str0 += r'$\theta$ = %.2f$^{\circ}$' % (popt[5]*180.0/np.pi)
    else:
        str0 = 'FWHM: saturated'

    t_ax.annotate(str0, [0.05,0.05], xycoords='axes fraction',
                  ha='left', va='bottom', fontsize=12)
    # print '## FWHM : ', FWHMx, FWHMy
#enddef

def get_slit_trace(infile): #, xmin, xmax):
    # Mod on 04/04/2017, aesthetics, fix bug
    # Mod on 04/04/2017, handle CRs affecting trace
    # Mod on 06/04/2017. Slit trace can be lost in median. Using average

    im0, hdr0 = fits.getdata(infile, header=True)

    # + on 04/04/2017
    im0_clean = cosmicray_lacosmic(im0, sigclip=10)
    im0_clean = im0_clean[0]

    # Bug: Mod on 06/04/2017
    y_avg0 = np.average(im0_clean, axis=1) # Mod on 04/04/2017
    cen0   = (np.where(y_avg0 == np.max(y_avg0))[0])[0]

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

    # + on 05/04/2017
    i_y_grp = group_index(i_y, find_max=True)
    i_x_grp = group_index(i_x, find_max=True)

    # + on 01/04/2017, Mod on 05/04/2017
    y_min, y_max = np.min(i_y_grp), np.max(i_y_grp)
    x_min, x_max = np.min(i_x_grp), np.max(i_x_grp)
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
    # print xcen, ycen

    xcen += pos[0]-find_size2d[1].value/2.0
    ycen += pos[1]-find_size2d[0].value/2.0
    # print xcen, ycen
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
    Modified by Chun Ly, 04-05 April 2017
     - Adjust greyscale limits to handle slit image (make it black),
       and faint sources
     Use find_gnirs_window_mean to find center
    Modified by Chun Ly, 05 April 2017
     - Handle alignment sequences with more than just 4 frames
     - Handle excess subplots for individual PDF pages (remove axes)
     - Compute seeing FWHM for acquisition images
    Modified by Chun Ly, 06 April 2017
     - Get coordinates for slit in cutout
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
            nrows = 2 # np.ceil(len(t_idx)/ncols)
            ncols, nrows = np.int(ncols), np.int(nrows)

            # Mod on 05/04/2017
            if len(t_idx) <= nrows * ncols:
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

                # + on 05/04/2017
                if len(t_idx) > (nrows*ncols):
                    if jj % (nrows * ncols) == 0:
                        fig, ax_arr = plt.subplots(nrows=nrows, ncols=ncols)

                im0  = fits.getdata(t_files[jj])
                hdr0 = fits.getheader(t_files[jj], ext=0) # Get WCS header

                # + 01/04/2017
                im0_clean = cosmicray_lacosmic(im0, sigclip=10)[0]

                cutout = Cutout2D(im0_clean, pos_cen, size2d,
                                  mode='partial', fill_value=np.nan)

                t_col, t_row = jj % ncols, (jj / ncols) % nrows

                # Mod on 04/04/2017 to handle bright and faint stars
                max0 = np.max(cutout.data)

                # Compute median within GNIRS window
                # + on 04-05/04/2017
                temp       = im0_clean[-50:-1,:]
                bgd0, sig0 = np.median(temp), np.std(temp)
                idx_y, idx_x = np.where(im0_clean > (bgd0 + 5*sig0))
                med0 = np.median(im0_clean[idx_y,idx_x])
                print '## max0 : ', max0, ' med0 : ', med0
                if max0 > 50000:
                    z1, z2 = zscale.get_limits(cutout.data)
                    z2 = max0 # Change for better stretch for telluric star
                else:
                    if ('Acq_' not in tab0['slit'][jj_idx]) and \
                       (tab0['exptime'][jj_idx] == 3):
                        # First frame that will show the longslit
                        z1, z2 = 0.0, 0.5*max0
                    else:
                        # This should handle faint and bright stars
                        z1, z2 = 0.5*med0, max0

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
                axins = zoomed_inset_axes(t_ax, 6, loc=4)
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

                # Compute FWHM of alignment star | + on 05/04/2017
                if ('Acq_' not in tab0['slit'][jj_idx]) and \
                   (tab0['exptime'][jj_idx] == 3):
                    log.info('## No source in slit : '+tab0['filename'][jj_idx])
                else:
                    # + on 06/04/2017
                    c_size2d     = u.Quantity((40, 40), u.pixel)
                    c_slit_x0    = slit_x0 - (c_xcen-c_size2d[1].value/2.0)
                    c_slit_y0_lo = slit_y0_lo - (c_ycen-c_size2d[0].value/2.0)
                    c_slit_y0_hi = slit_y0_hi - (c_ycen-c_size2d[0].value/2.0)
                    im0_crop = Cutout2D(cutout.data, (c_xcen,c_ycen), c_size2d,
                                        mode='partial', fill_value=np.nan)
                    gauss2d_fit(im0_crop.data, hdr0, t_ax, c_slit_x0,
                                c_slit_y0_lo, c_slit_y0_hi) # Mod on 06/04/2017

                # Write each page separately | + on 05/04/2017
                if len(t_idx) > (nrows*ncols):
                    # Mod later on 05/04/2017 to handle excess subplots
                    if jj == len(t_idx)-1:
                        rem0 = len(t_idx) % (nrows*ncols) # remainder
                        if rem0 != 0:
                            for rr in range(rem0,nrows*ncols,1):
                                t_col, t_row = rr % ncols, (rr / ncols) % nrows
                                ax_arr[t_row,t_col].axis('off')

                    if (jj % (nrows * ncols) == nrows*ncols-1) or \
                       (jj == len(t_idx)-1):
                        subplots_adjust(left=0.02, bottom=0.02, top=0.95,
                                        right=0.98, wspace=0.02, hspace=0.02)
                        fig.set_size_inches(11,8)
                        fig.savefig(pp, format='pdf')
            #endfor

            # Mod on 05/04/2017
            if len(t_idx) <= nrows * ncols:
                # Mod later on 05/04/2017 to handle excess subplots
                for rr in range(len(t_idx),nrows*ncols):
                    t_col, t_row = rr % ncols, (rr / ncols) % nrows
                    ax_arr[t_row,t_col].axis('off')

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
