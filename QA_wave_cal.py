"""
QA_wave_cal
===========

A set of code to provide plots illustrating wavelength calibration solutions
"""

import sys, os

from chun_codes import systime

from os.path import exists
import glob

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import fits
from astropy.io import ascii as asc # + on 16/06/2017
from astropy import log

# + on 19/05/2017
from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()
from astropy.visualization.mpl_normalize import ImageNormalize

import aplpy # + on 19/05/2017

import glog # + on 09/01/2018

from pylab import subplots_adjust
bbox_props = dict(boxstyle="square,pad=0.15", fc="w", alpha=0.5, ec="none")

co_dirname = os.path.dirname(__file__)

def arc_check(path, arcs=[''], out_pdf='', stack=False, silent=False, verbose=True):

    '''
    Generate plot illustrating wavelength calibration from IRAF database

    Parameters
    ----------
    path : str
      Full path to where output PDF and FITS file are located. Must end
      with a '/'

    arcs : str or list (Optional)
      List of raw filenames for the arc data (e.g., 'N20170101S0111.fits')

    out_pdf : str
      Filename for output PDF. Do NOT include full path

    stack : boolean
      Indicate whether to use stacked arc data. Default: False

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 18 May 2017
    Modified by Chun Ly, 19 May 2017
     - Fix bug with scatter plot (wrong indexing)
     - Fix bug with l_val0 and l_val (those were flipped)
     - Draw dashed lines connecting the points
     - Adjust margins, add labels

     - Adjust image scale (Using zscale)
     - Annotate each page with the dataframe information
     - Label lines
     - Aesthetic changes for plotting

     - Aesthetic changes for plotting
       - subplots_adjust
       - Better labeling of lines
       - xlim, ylim set
       - Adjust PDF paper size

    Modified by Chun Ly, 11 November 2017
     - Added stack keyword option to use stacked arc products

    Modified by Chun Ly, 15 November 2017
     - Bug fix: Incorrect replace

    Modified by Chun Ly, 16 November 2017
     - Change prefix: rnc to rbnc

    Modified by Chun Ly,  9 January 2018
     - Import glog and call for stdout and ASCII logging
    '''

    # + on 09/01/2018
    logfile  = path0+'QA_wave_cal.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin arc_check : '+systime())

    # Mod on 11/11/2017
    if stack == True:
        arcs      = ['arc_stack.fits'] # + on 15/11/2017
        rnc_files = [path + 'arc_stack.fits']
        d_files = [path+'database/idwarc_stack_SCI_1_']
    else:
        if arcs[0] == '':
            arc_list = path + 'arc.lis'
            if silent == False: mylogger.info('Reading : '+arc_list)
            arcs = np.loadtxt(arc_list, dtype=type(str))


        # Mod on 16/11/2017
        rnc_files = [path+'rbnc'+file0 for file0 in arcs]
        d_files   = [path+'database/idwrbnc'+file0.replace('.fits','_SCI_1_') \
                     for file0 in arcs]

    n_arcs = len(rnc_files) # Mod on 11/11/2017

    out_pdf = path+'arc_check.pdf' if out_pdf == '' else path+out_pdf
    pp = PdfPages(out_pdf)

    for nn in xrange(n_arcs):
        hdu0 = fits.open(rnc_files[nn])
        im0  = hdu0['SCI'].data

        fig, ax = plt.subplots()
        z1, z2 = zscale.get_limits(im0)
        norm = ImageNormalize(vmin=z1, vmax=z2)
        ax.imshow(im0, cmap='Greys', origin='lower', norm=norm)

        f0   = open(d_files[nn], 'r')
        str0 = f0.readlines()

        beg_mark  = [ii for ii in xrange(len(str0)) if 'begin' in str0[ii]]
        feat_mark = [ii for ii in xrange(len(str0)) if 'features' in str0[ii]]
        func_mark = [ii for ii in xrange(len(str0)) if 'function' in str0[ii]]

        n_cols     = len(beg_mark)
        n_features = [np.int(val.split('\t')[2].replace('\n','')) for val in
                      np.array(str0)[feat_mark]]
        n_max      = max(n_features)
        x_cols, y_val = np.zeros(n_cols), np.zeros((n_cols,n_max))
        l_val0, l_val = np.zeros((n_cols,n_max)), np.zeros((n_cols,n_max))
        for cc in xrange(n_cols):
            temp = str0[beg_mark[cc]].split(' ')[1]
            # Mod on 16/11/2017
            temp = temp.replace('wrbnc'+arcs[nn].replace('.fits','[SCI,1]['),'')
            temp = temp.replace('w'+arcs[nn].replace('.fits','[SCI,1]['),'') # + on 15/11/2017
            x_cols[cc] = np.int(temp.split(',')[0])

            for rr in xrange(n_features[cc]): #beg_mark[cc]+1, func_mark[cc]):
                temp1 = str0[feat_mark[cc]+1+rr].replace('       ','')
                temp1 = temp1.split(' ')
                temp1 = [val for val in temp1 if val != '']
                y_val[cc,rr]  = np.float(temp1[1])
                try:
                    l_val[cc,rr]  = np.float(temp1[2])
                except ValueError:
                    pass

                try:
                    l_val0[cc,rr] = np.float(temp1[3])
                except ValueError:
                    pass
            #endfor
            x_temp = np.repeat(x_cols[cc],n_features[cc])
            y_temp = y_val[cc,0:n_features[cc]]
            ax.scatter(x_temp, y_temp, facecolor='none', edgecolor='red',
                       marker='o', s=25, alpha=0.5)
        #endfor

        line_list = list(set(l_val0.reshape(n_cols*n_max).tolist()))
        for ll in range(len(line_list)):
            if line_list[ll] != 0:
                x_idx, y_idx = np.where(l_val0 == line_list[ll])
                sort0 = np.argsort(x_cols[x_idx])
                x1 = x_cols[x_idx][sort0]
                y1 = y_val[x_idx,y_idx][sort0]
                ax.plot(x1, y1, 'r--', alpha=0.5)
                ax.annotate(str(line_list[ll]), [690,np.average(y1)],
                            ha='right', va='center', fontsize=10,
                            color='red', bbox=bbox_props)
        ax.set_xlabel('X [pixels]')
        ax.set_ylabel('Y [pixels]')
        ax.set_xlim([0,690])
        ax.set_ylim([0,1025])
        ax.minorticks_on()
        subplots_adjust(left=0.095, bottom=0.025, top=0.995, right=0.99)
        str0 = (rnc_files[nn]+'\n'+d_files[nn]).replace(path,'')
        ax.annotate(str0, [0.025,0.975], xycoords='axes fraction', ha='left',
                    va='top', fontsize=14, bbox=bbox_props)
        ax.set_title(path)
        fig.set_size_inches(7.65,10.5) # Scale proportionally
        fig.savefig(pp, format='pdf')
    #endfor

    pp.close()
    if silent == False: mylogger.info('### End arc_check : '+systime())
#enddef

def OH_check(path, objs='', out_pdf='', skysub=False, silent=False,
             verbose=True):

    '''
    Generate plot illustrating expected location of OH skyline to
    check wavelength calibration

    Parameters
    ----------

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    skysub : boolean
      Display skysubtracted or un-skysubtracted images. Default: False

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 19 May 2017
    Modified by Chun Ly, 20 May 2017
     - Draw OH night skyline emission on plots
    Modified by Chun Ly, 2 June 2017
     - Add skysub keyword option to use either sky-subtracted
       images or those without skysubtraction (OH skylines more visible)
    Modified by Chun Ly, 26 June 2017
     - Plot more OH skylines (change threshold for OH_mark)
    Modified by Chun Ly, 16 November 2017
     - Change prefix: tfrnc to tfrbnc
    Modified by Chun Ly,  9 January 2018
     - Import glog and call for stdout and ASCII logging
    '''

    # + on 09/01/2018
    logfile  = path0+'QA_wave_cal.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin OH_check : '+systime())

    # + on 20/05/2017
    OH_file = co_dirname+'/rousselot2000.dat'
    if exists(OH_file):
        if silent == False: mylogger.info('Reading : '+OH_file)
        OH_data  = np.loadtxt(OH_file)
        OH_lines = OH_data[:,0]
        OH_int   = OH_data[:,1]

    # Mod on 02/06/2017
    if objs == '':
        obj_list = path + 'obj.lis' if skysub == True else \
                   path + 'obj.OH.lis'

        if silent == False: mylogger.info('Reading : '+obj_list)
        objs = np.loadtxt(obj_list, dtype=type(str))

    n_obj = len(objs)

    tfrnc_files = [path+'tfrbnc'+file0 for file0 in objs] # Mod on 16/11/2017

    # + on 20/05/2017
    chk = [file0 for file0 in tfrnc_files if exists(file0) == True]
    if len(chk) == 0:
        mylogger.warn('Files not found!!!')
        mylogger.warn(', '.join(file0))
        mylogger.warn('Exiting!!!')
        return

    # Mod on 02/06/2017
    if out_pdf == '':
        out_pdf = path+'OH_check.pdf' if skysub == True else \
                  path+'OH_check.raw.pdf'
    else:
        out_pdf = path+out_pdf

    pp = PdfPages(out_pdf)

    for nn in xrange(n_obj):
        if exists(tfrnc_files[nn]): # Mod on 20/05/2017
            hdu0 = fits.open(tfrnc_files[nn])
            im0  = hdu0['SCI'].data

            gc0 = aplpy.FITSFigure(hdu0, hdu='SCI', figsize=(7.65,10.5))

            z1, z2 = zscale.get_limits(im0)
            gc0.show_grayscale(invert=True, vmin=z1, vmax=z2)
            gc0.set_tick_color('black')

            # + on 20/05/2017
            hdr0     = hdu0['SCI'].header
            crval2   = hdr0['CRVAL2']
            cdelt2   = hdr0['CDELT2']
            npix     = hdr0['NAXIS2']
            lamb_max = crval2 + cdelt2*npix
            OH_mark  = np.where((OH_lines >= crval2) & (OH_lines <= lamb_max))[0]
            # Mod on 26/06/2017
            OH_mark  = np.where(OH_int >= 0.005*np.max(OH_int[OH_mark]))[0]
            line_list = []
            for ll in OH_mark:
                line_list.append(np.array([[0,700], [OH_lines[ll],OH_lines[ll]]]))
            gc0.show_lines(line_list, color='red', alpha=0.5, linewidth=1.0,
                           linestyle='dashed')

            #subplots_adjust(left=0.095, bottom=0.025, top=0.995, right=0.99)
            str0 = tfrnc_files[nn].replace(path,'')
            gc0.add_label(0.025, 0.975, str0, color='red', relative=True,
                          ha='left', va='top', weight='medium', fontsize=14,
                          bbox=bbox_props)
            gc0.set_axis_labels(xlabel='X [pixels]', ylabel=r'Wavelength ($\AA$)')
            #gc0.set_axis_label_rotation(90)
            gc0.savefig(pp, format='pdf')
        #endif
    #endfor

    if silent == False: mylogger.info('Writing : '+out_pdf)
    pp.close()
    if silent == False: mylogger.info('### End OH_check : '+systime())
#enddef

def arc_check2(path, arcs=[''], out_pdf='', stack=False, silent=False, verbose=True):

    '''
    Generate plot illustrating expected location of arc lines to
    check wavelength calibration

    Parameters
    ----------
    path : str
      Full path to where output PDF and FITS file are located. Must end
      with a '/'

    arcs : str or list (Optional)
      List of raw filenames for the arc data (e.g., 'N20170101S0111.fits')

    stack : boolean
      Indicate whether to include stacked arc data. Default: False

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 16 June 2017
     - Started as a copy of OH_check
    Modified by Chun Ly, 17 July 2017
     - Fix bug with arcs handling
    Modified by Chun Ly, 12 November 2017
     - Added stack keyword option to include stacked arc products
    Modified by Chun Ly, 16 November 2017
     - Change prefix: tfrnc to tfrbnc
    '''

    if silent == False: log.info('### Begin arc_check2 : '+systime())

    arc_file = co_dirname+'/argon.dat'
    if exists(arc_file):
        if silent == False: log.info('### Reading : '+arc_file)
        arc_data   = asc.read(arc_file, format='commented_header')
        arc_lines  = arc_data['Line'].data
        arc_source = arc_data['Source'].data

    if arcs[0] == '':
        arcs_list = path + 'arc.lis'

        if silent == False: log.info('### Reading : '+arcs_list)
        arcs = np.loadtxt(arcs_list, dtype=type(str))

    tfrnc_files = [path+'tfrbnc'+file0 for file0 in arcs] # Mod on 16/11/2017

    # + on 12/11/2017
    if stack == True:
        arcs0        = [path + 'arc_stack.fits']
        tfrnc_files0 = [path + 'tfarc_stack.fits']

        arcs         = arcs0 + arcs
        tfrnc_files  = tfrnc_files0 + tfrnc_files

    chk = [file0 for file0 in tfrnc_files if exists(file0) == True]
    if len(chk) == 0:
        log.warn('### Files not found!!!')
        log.warn('### '+', '.join(file0))
        log.warn('### Exiting!!!')
        return

    n_arcs = len(arcs) # Moved lower on 12/11/2017

    out_pdf = path+'arc_check2.pdf' if out_pdf == '' else path+out_pdf
    pp = PdfPages(out_pdf)

    for nn in xrange(n_arcs):
        if exists(tfrnc_files[nn]): # Mod on 20/05/2017
            hdu0 = fits.open(tfrnc_files[nn])
            im0  = hdu0['SCI'].data

            gc0 = aplpy.FITSFigure(hdu0, hdu='SCI', figsize=(7.65,10.5))

            z1, z2 = zscale.get_limits(im0)
            gc0.show_grayscale(invert=True, vmin=z1, vmax=z2)
            gc0.set_tick_color('black')

            # + on 20/05/2017
            hdr0     = hdu0['SCI'].header
            crval2   = hdr0['CRVAL2']
            cdelt2   = hdr0['CDELT2']
            npix     = hdr0['NAXIS2']
            lamb_max = crval2 + cdelt2*npix
            arc_mark = np.where((arc_lines >= crval2) & (arc_lines <= lamb_max))[0]
            line_list = []
            for ll in arc_mark:
                line_list.append(np.array([[0,700], [arc_lines[ll],arc_lines[ll]]]))
            gc0.show_lines(line_list, color='red', alpha=0.5, linewidth=1.0,
                           linestyle='dashed')

            str0 = tfrnc_files[nn].replace(path,'')
            gc0.add_label(0.025, 0.975, str0, color='red', relative=True,
                          ha='left', va='top', weight='medium', fontsize=14,
                          bbox=bbox_props)
            gc0.set_axis_labels(xlabel='X [pixels]', ylabel=r'Wavelength ($\AA$)')
            gc0.savefig(pp, format='pdf')
        #endif
    #endfor

    if silent == False: log.info('### Writing : '+out_pdf)
    pp.close()
    if silent == False: log.info('### End arc_check2 : '+systime())
#enddef
