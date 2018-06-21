"""
QA_wave_cal
===========

A set of code to provide plots illustrating wavelength calibration solutions
"""

import sys, os

from pyraf import iraf

from chun_codes import systime

from os.path import exists
import glob

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from scipy.optimize import curve_fit # + on 05/06/2018
from IQ_plot import gauss1d
from get_OH_centers import gauss_multi, n_multi

from astropy.io import fits
from astropy.io import ascii as asc # + on 16/06/2017
from astropy import log
from astropy.table import Table

# + on 19/05/2017
from astropy.visualization import ZScaleInterval
zscale = ZScaleInterval()
from astropy.visualization.mpl_normalize import ImageNormalize

import aplpy # + on 19/05/2017

import glog # + on 09/01/2018

from pylab import subplots_adjust
bbox_props = dict(boxstyle="square,pad=0.15", fc="w", alpha=0.5, ec="none")

iraf.gemini(_doprint=0)
iraf.gemini.gnirs(_doprint=0)

log.info("Unlearning tasks")
iraf.gemini.unlearn()
iraf.gemini.gemtools.unlearn()
iraf.gemini.gnirs.unlearn()

iraf.set(stdimage="imt4096")
# iraf.gemini.nsheaders("gnirs")

co_dirname = os.path.dirname(__file__)

xorder = 3 # nsfitcoords fitting order along x

def get_database_model(path, source, silent=False, verbose=True):
    '''
    Determine fitting function and order to pass to nsfitcoords

    Parameters
    ----------
    path : str
      Full path to where output PDF and FITS file are located. Must end
      with a '/'

    source : str
      Either 'arc' or 'OH'

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True
    '''

    logfile  = path+'QA_wave_cal.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin get_database_model : '+systime())

    ref_file = source+'_stack.fits'

    if source == 'arc': dbase = 'database'
    if source == 'OH':  dbase = 'database_OH'

    d_file  = '%s%s/idw%s' % (path, dbase, ref_file.replace('.fits','_SCI_1_'))

    mylogger.info("Reading : "+d_file)

    f0   = open(d_file, 'r')
    str0 = f0.readlines()

    ii_beg  = [ii for ii in xrange(len(str0)) if 'begin' in str0[ii]]
    ii_feat = [ii for ii in xrange(len(str0)) if 'features' in str0[ii]]
    ii_func = [ii for ii in xrange(len(str0)) if 'function' in str0[ii]]
    ii_ord  = [ii for ii in xrange(len(str0)) if 'order' in str0[ii]]

    func0  = str0[ii_func[0]].split(' ')[-1].replace('\n','')
    order0 = np.int(str0[ii_ord[0]].split(' ')[-1].replace('\n',''))

    mylogger.info('Function for '+source+' line fitting : '+func0)
    mylogger.info('Order of '+source+' line fitting : %i' % order0)
    f0.close()

    if silent == False: mylogger.info('### End get_database_model : '+systime())

    return func0, order0
#enddef

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
    logfile  = path+'QA_wave_cal.log'
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
             verbose=True, cross_check=False):

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

    cross_check : boolean
      Check arc-based wavelength calibration against OH skylines. Default: False

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
    Modified by Chun Ly, 31 May 2018
     - Add cross_check keyword; Update input/outputs for cross_check == True
     - Bug fix: n_obj for cross_check == True
    '''

    # + on 09/01/2018
    logfile  = path+'QA_wave_cal.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin OH_check : '+systime())

    # + on 20/05/2017
    OH_file = co_dirname+'/rousselot2000.dat'
    if exists(OH_file):
        if silent == False: mylogger.info('Reading : '+OH_file)
        OH_data  = np.loadtxt(OH_file)
        OH_lines = OH_data[:,0]
        OH_int   = OH_data[:,1]

    # Mod on 02/06/2017, 31/05/2018
    if cross_check == False:
        if objs == '':
            obj_list = path + 'obj.lis' if skysub == True else \
                       path + 'obj.OH.lis'

            if silent == False: mylogger.info('Reading : '+obj_list)
            objs = np.loadtxt(obj_list, dtype=type(str))

        n_obj = len(objs)

        tfrnc_files = [path+'tfrbnc'+file0 for file0 in objs] # Mod on 16/11/2017
    else:
        tfrnc_files = [path+'tfOH_stack_arc.fits']
        n_obj = len(tfrnc_files)

    # + on 20/05/2017
    chk = [file0 for file0 in tfrnc_files if exists(file0) == True]
    if len(chk) == 0:
        mylogger.warn('Files not found!!!')
        mylogger.warn(', '.join(file0))
        mylogger.warn('Exiting!!!')
        return

    # Mod on 02/06/2017, 31/05/2018
    if cross_check == False:
        if out_pdf == '':
            out_pdf = path+'OH_check.pdf' if skysub == True else \
                      path+'OH_check.raw.pdf'
        else:
            out_pdf = path+out_pdf
    else:
        if out_pdf == '':
            out_pdf = path+'OH_check_arc.pdf'
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

def arc_check2(path, arcs=[''], out_pdf='', stack=False, cross_check=False,
               silent=False, verbose=True):

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

    cross_check : boolean
      Check OH-based wavelength calibration against arc lines. Default: False

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
    Modified by Chun Ly,  9 January 2018
     - Import glog and call for stdout and ASCII logging
    Modified by Chun Ly, 31 May 2018
     - Add cross_check keyword; Update input/outputs for cross_check == True
    '''

    # + on 09/01/2018
    logfile  = path+'QA_wave_cal.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin arc_check2 : '+systime())

    arc_file = co_dirname+'/argon.dat'
    if exists(arc_file):
        if silent == False: mylogger.info('Reading : '+arc_file)
        arc_data   = asc.read(arc_file, format='commented_header')
        arc_lines  = arc_data['Line'].data
        arc_source = arc_data['Source'].data

    if arcs[0] == '':
        arcs_list = path + 'arc.lis'

        if silent == False: mylogger.info('Reading : '+arcs_list)
        arcs = np.loadtxt(arcs_list, dtype=type(str))

    tfrnc_files = [path+'tfrbnc'+file0 for file0 in arcs] # Mod on 16/11/2017

    # + on 12/11/2017
    if stack == True:
        arcs0        = [path + 'arc_stack.fits']
        tfrnc_files0 = [path + 'tfarc_stack.fits']

        arcs         = arcs0 + arcs
        tfrnc_files  = tfrnc_files0 + tfrnc_files

    # + on 31/05/2018
    if cross_check == True:
        arcs         = [path + 'arc_stack_OH.fits']
        tfrnc_files  = [path + 'tfarc_stack_OH.fits']

    chk = [file0 for file0 in tfrnc_files if exists(file0) == True]
    if len(chk) == 0:
        mylogger.warn('Files not found!!!')
        mylogger.warn(', '.join(file0))
        mylogger.warn('Exiting!!!')
        return

    n_arcs = len(arcs) # Moved lower on 12/11/2017

    # Mod on 31/05/2018
    if cross_check == False:
        out_pdf = path+'arc_check2.pdf' if out_pdf == '' else path+out_pdf
    else:
        out_pdf = path+'arc_check_OH.pdf' if out_pdf == '' else path+out_pdf

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

    if silent == False: mylogger.info('Writing : '+out_pdf)
    pp.close()
    if silent == False: mylogger.info('### End arc_check2 : '+systime())
#enddef

def residual_wave_cal(path, dataset='', cal='', silent=False, verbose=True):
    '''
    Plot residuals of arc/OH lines against OH/arc dataset

    Parameters
    ----------
    path : str
      Full path to where output PDF and FITS file are located. Must end
      with a '/'

    dataset : str
      Either: 'arc' or 'OH'

    cal : str
      Either 'arc' or 'OH

    To use OH calib on OH lines:   dataset='OH',  cal='OH'
    To use OH calib on arc lines:  dataset='arc', cal='OH'
    To use arc calib on arc lines: dataset='arc', cal='arc'
    To use arc calib on OH lines:  dataset='OH',  cal='arc'

    silent : boolean
      Turns off stdout messages. Default: False

    verbose : boolean
      Turns on additional stdout messages. Default: True

    Returns
    -------

    Notes
    -----
    Created by Chun Ly, 4 June 2018
    Modified by Chun Ly, 5 June 2018
     - Define and read in OH/arc line database
     - Import curve_fit
     - Import gauss1d from IQ_plot
     - Call curve_fit to fit Gaussian profiles to line
    Modified by Chun Ly, 7 June 2018
     - Plot residuals of wavelength solution
     - Set minimum line peak brightness to fit
     - Limit plotting to non-zero values for central wavelength
     - Plot rms for each line
     - Write average array to FITS file
     - Plot aesthetics
     - Move mylogger.info() call into else statement
     - Handle dataset == cal cases
     - Save average and rms for each arc/OH line
     - Handle dataset == cal cases for output PDF file
     - Plot aesthetics: Draw vertical lines for location of lines
     - Handle dataset == cal cases for output FITS file
     - If statement for lines with fits
     - Write ASCII with averages and RMS file
    Modified by Chun Ly, 19 June 2018
     - Determine and plot average, median, and rms for each column
     - Multi-page PDF: Separate plot illustrating difference from model
       along longslit
     - Plot sigma on second page, Exclude no value cases
     - Add ax.annotation in upper left of second page
     - Call get_database_model, include fitting model
     - Overwrite ASCII tables in asc.write
    Modified by Chun Ly, 20 June 2018
     - Use FITS file if it exists, otherwise generate
     - Use sigma_clipped_stats to remove outliers in statistics
     - Remove outlier based on offset amount, not using sigma clipping
     - Tab issue: Move fits.writeto out of for loop
     - Compute stats for histogram distributions
     - try/except for OH skyline failure (RunTimeError)
    Modified by Chun Ly, 21 June 2018
     - Use convolved Rousselot2000 table instead of original
     - Read in OH npz file for line grouping
     - Define use_lines array for arc case
     - Handle multi-line fitting (works with dataset='arc')
     - Handle multi-line fitting for dataset='OH'
     - Include lower and upper bounds, Move multi-line fitting later in code
     - Normalize spectrum to peak, Do not set bounds for curve_fit()
    '''

    logfile  = path+'QA_wave_cal.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin residual_wave_cal : '+systime())

    if dataset != cal:
        infile = path+'tf'+dataset+'_stack_'+cal+'.fits'
    else:
        infile = path+'tf'+dataset+'_stack.fits'

    if not exists(infile):
        mylogger.warn('File not found : '+infile)
        mylogger.warn('Exiting!!!')
        return
    else:
        if silent == False: mylogger.info('Reading : '+infile)
        cal_2D, cal_hdr = fits.getdata(infile, extname='SCI', header=True)

        NX = cal_2D.shape[1] # spatial direction
        NY = cal_2D.shape[0] # dispersion direction

        lam0 = cal_hdr['CRVAL2']
        dlam = cal_hdr['CD2_2']
        wave0 = lam0 + dlam * np.arange(NY)

        bins_mid = np.arange(13,NX,10)
        n_bins = len(bins_mid)

        if dataset != cal:
            out_fits = path+'wave_cal_resid_'+dataset+'_'+cal+'.fits'
        else:
            out_fits = path+'wave_cal_resid_'+dataset+'.fits'

        if not exists(out_fits):
            avg_arr = np.zeros( (NY, n_bins) )
            for ii in range(n_bins):
                start = np.max([0,bins_mid[ii]-1 - 5])
                stop  = np.min([NX-1,bins_mid[ii]-1 + 5])
                avg_arr[:,ii] = np.average(cal_2D[:,start:stop], axis=1)

            if silent == False: mylogger.info('Writing : '+out_fits)
            fits.writeto(out_fits, avg_arr, overwrite=True)
        else:
            if silent == False: mylogger.info('Reading : '+out_fits)
            avg_arr = fits.getdata(out_fits)

        if dataset == 'OH':
            cal_ref_file = path+'rousselot2000_convl.dat'

        if dataset == 'arc':
            cal_ref_file = co_dirname+'/argon.dat'

        if silent == False: mylogger.info('Reading : '+cal_ref_file)
        cal_line_data = asc.read(cal_ref_file, format='no_header')
        cal_lines     = cal_line_data['col1'].data

        skip = np.zeros(len(cal_lines))

        if dataset == 'OH':
            npz_file = path+'rousselot2000_convl.npz'
            #mylogger.info('Reading : '+npz_file)
            npz_tab = np.load(npz_file)
            use_lines = npz_tab['use_lines']
            for gg in range(len(use_lines)):
                matches = [xx for xx in range(len(cal_lines)) if
                           cal_lines[xx] in use_lines[gg]]
                if len(matches) > 1:
                    skip[np.array(matches[1:])] = 1

        in_spec = np.where((cal_lines >= wave0[0]) &
                           (cal_lines <= wave0[-1]))[0]
        cal_line_data = cal_line_data[in_spec]
        cal_lines     = cal_line_data['col1'].data
        skip          = skip[in_spec]

        if dataset == 'OH':
            cal_lines_int = cal_line_data['col2'].data

            # Remove lines from use_lines if outside coverage | + on 21/06/2018
            n_in_spec = np.zeros(len(use_lines))
            for gg in range(len(use_lines)):
                in_spec = np.where((use_lines[gg] >= wave0[0]) &
                                   (use_lines[gg] <= wave0[-1]))[0]
                n_in_spec[gg] = len(in_spec)
                if len(in_spec) != len(use_lines[gg]) and len(in_spec) > 0:
                    use_lines[gg] = np.array(use_lines[gg])[in_spec].tolist()
            ignore = np.where(n_in_spec == 0)[0]
            if len(ignore) > 0:
                use_lines = np.delete(use_lines, ignore)

        n_lines = len(cal_lines)

        if dataset == 'arc':
            use_lines = []
            for ll in range(n_lines): use_lines.append([cal_lines[ll]])
            use_lines = np.array(use_lines)

        cen_arr = np.zeros( (n_lines,n_bins) )

        if dataset != cal:
            pdf_file = path+'wave_cal_resid_'+dataset+'_'+cal+'.pdf'
        else:
            pdf_file = path+'wave_cal_resid_'+dataset+'.pdf'

        fig, ax = plt.subplots()

        diff_avg = np.zeros(n_lines)
        diff_rms = np.zeros(n_lines)
        diff_num = np.zeros(n_lines, dtype=np.int)

        u_l_ii = 0
        for ll in range(n_lines):
            mylogger.info('ll=%i, u_l_ii=%i' % (ll, u_l_ii))
            ax.axvline(cal_lines[ll], color='red', linestyle='dashed',
                       linewidth=0.25, zorder=1)

            if not skip[ll]:
                w_min = min(use_lines[u_l_ii])-10
                w_max = max(use_lines[u_l_ii])+10
                z_idx = np.where((wave0 >= w_min) & (wave0 <= w_max))[0]
                x0 = wave0[z_idx]

                for ii in range(n_bins):
                    y0 = avg_arr[z_idx,ii]
                    y0 /= max(y0)
                    if len(use_lines[u_l_ii]) == 1:
                        p0 = [0.0, 1.0, cal_lines[ll], 2.0]
                        if p0[1] > 10:
                            try:
                                popt, pcov = curve_fit(gauss1d, x0, y0, p0=p0)
                                cen_arr[ll,ii] = popt[2]
                            except RuntimeError:
                                print ll, ii, p0
                    else:
                        # Set-up initial curve_fit guess for multi lines
                        if len(use_lines[u_l_ii]) > 1:
                            matches = np.array([xx for xx in range(n_lines) if \
                                                cal_lines[xx] in use_lines[u_l_ii]])
                            ratio0  = cal_lines_int[matches]/max(cal_lines_int[matches])
                            t_peak0 = ratio0.tolist()
                            t_lines = list(use_lines[u_l_ii])
                            t_sig   = [2.0] * len(t_lines)

                            t_peak0 += np.zeros(n_multi-len(t_lines)).tolist()
                            t_sig   += np.zeros(n_multi-len(t_lines)).tolist()
                            t_lines += np.zeros(n_multi-len(t_lines)).tolist()

                            p0 = t_peak0
                            p0 += t_lines
                            p0 += t_sig

                            p0 = np.array(p0)
                            low_bound = tuple([0] * n_multi) + \
                                        tuple(p0[n_multi:2*n_multi]-0.5) + \
                                        tuple([0] * n_multi)
                            up_bound  = tuple(p0[0:n_multi]*1.25+0.1) + \
                                        tuple(p0[n_multi:2*n_multi]+0.5) + \
                                        tuple(p0[2*n_multi:]+1)

                            try:
                                popt, pcov = curve_fit(gauss_multi, x0, y0, p0=p0)
                                #bounds=(low_bound, up_bound))
                                t_loc = popt[n_multi:2*n_multi]
                                cen_arr[matches,ii] = t_loc[np.where(t_loc != 0)[0]]
                            except RuntimeError:
                                pass #print 'fail : ', ll, ii, p0, popt

                #endfor
                u_l_ii += 1

                good = np.where(cen_arr[ll] != 0)[0]
                if len(good) > 0:
                    x_test = cen_arr[ll][good]
                    diff = x_test - cal_lines[ll]
                    ax.scatter(x_test, diff, marker='o', s=5, edgecolor='none',
                               facecolor='black', alpha=0.5, zorder=2)
                    diff_rms[ll], diff_avg[ll] = np.std(diff), np.average(diff)
                    diff_num[ll] = len(good)

            #endif
        #endfor
        ax.scatter(cal_lines, diff_avg, marker='o', s=40, facecolor='blue',
                   edgecolor='none', alpha=0.5, zorder=3)
        ax.errorbar(cal_lines, diff_avg, yerr=diff_rms, ecolor='blue',
                    capsize=2.0, elinewidth=2, fmt=None, alpha=0.5, zorder=4)

        ax.minorticks_on()
        ax.set_xlabel(r'Wavelengths [$\AA$]')
        ax.set_ylabel(r'Difference from reference values [$\AA$]')

        func0, order0 = get_database_model(path, cal)
        an_txt  = 'dataset: %s  calib: %s\n' % (dataset, cal)
        an_txt += 'func: %s  yorder: %i  xorder: %i' % (func0, order0, xorder)

        ax.annotate(an_txt, [0.05,0.95], xycoords='axes fraction',
                    ha='left', va='top')
        plt.subplots_adjust(left=0.1, right=0.99, bottom=0.1, top=0.99)
        fig.savefig(pdf_file)

        asc_file = pdf_file.replace('.pdf', '.tbl')
        tab0 = Table([cal_lines, diff_num, diff_avg, diff_rms],
                     names=('Line','N','Avg','RMS'))
        asc.write(tab0, asc_file, format='fixed_width_two_line', overwrite=True)

        # Determine average, median, and rms for each column | + on 19/06/2018
        avg0 = np.zeros(n_bins)
        N0   = np.zeros(n_bins)
        med0 = np.zeros(n_bins)
        rms0 = np.zeros(n_bins)

        for nn in range(n_bins):
            t_diff   = cen_arr[:,nn] - cal_lines
            good = np.where((cen_arr[:,nn] != 0) &
                            (np.absolute(t_diff) <= 1.5))[0]
            if len(good) > 0:
                N0[nn]   = len(good)
                avg0[nn] = np.average(t_diff[good])
                med0[nn] = np.median(t_diff[good])
                rms0[nn] = np.std(t_diff[good])

        pdf_file2 = pdf_file.replace('.pdf', '.stat.pdf')
        pp = PdfPages(pdf_file2)

        fig, ax = plt.subplots(ncols=3)
        good = np.where(rms0 != 0)
        N1, b1, p1 = ax[0].hist(avg0[good], bins=10, align='mid', color='b',
                                alpha=0.5, edgecolor='none',
                                histtype='stepfilled')
        N2, b2, p2 = ax[1].hist(med0[good], bins=10, align='mid', color='m',
                                alpha=0.5, edgecolor='none',
                                histtype='stepfilled')
        N3, b3, p3 = ax[2].hist(rms0[good], bins=10, align='mid', color='k',
                                alpha=0.5, edgecolor='none',
                                histtype='stepfilled')

        ax[0].set_ylabel('N')
        ax[0].set_xlabel(r'Average [$\AA$]')
        ax[1].set_xlabel(r'Median [$\AA$]')
        ax[2].set_xlabel(r'$\sigma$ [$\AA$]')

        stat_avg = [np.average(avg0[good]), np.average(med0[good]),
                    np.average(rms0[good])]
        stat_med = [np.median(avg0[good]), np.median(med0[good]),
                    np.median(rms0[good])]
        stat_rms = [np.std(avg0[good]), np.std(med0[good]), np.std(rms0[good])]

        y_max = np.max([np.max(N1),np.max(N2),np.max(N3)])*1.1
        for aa in range(3):
            ax[aa].set_ylim([0,y_max])

            ax[aa].axvline(stat_avg[aa], linestyle='solid')
            ax[aa].axvline(stat_med[aa], linestyle='dashed')

            txt0 = ''
            str0 = [r'$<x>$', r'$\tilde{x}$', r'$\sigma$']
            val  = [stat_avg[aa], stat_med[aa], stat_rms[aa]]
            for ss,vv in zip(str0,val): txt0 += ss+(' = %.3f' % vv) +'\n'
            ax[aa].annotate(txt0, [0.95,0.95], xycoords='axes fraction',
                            ha='right', va='top')
        ax[1].set_yticklabels([])
        ax[2].set_yticklabels([])

        ax[1].tick_params(axis='y', direction='in')
        ax[2].tick_params(axis='y', direction='in')

        plt.subplots_adjust(left=0.08, right=0.99, bottom=0.12, top=0.96,
                            wspace=0.03)
        fig.set_size_inches(8,4)
        fig.savefig(pp, format='pdf')

        asc_file2 = pdf_file2.replace('.pdf', '.tbl')
        tab2 = Table([bins_mid, N0, avg0, med0, rms0],
                     names=('Column','N','Avg','Med','RMS'))
        tab2.write(asc_file2, format='ascii.fixed_width_two_line',
                   overwrite=True)

        fig, ax = plt.subplots()
        ax.scatter(bins_mid[good], avg0[good], marker='o', color='b',
                   alpha=0.5, edgecolor='none', label='Average')
        ax.scatter(bins_mid[good], med0[good], marker='o', color='m',
                   alpha=0.5, edgecolor='none', label='Median')
        ax.scatter(bins_mid[good], rms0[good], marker='o', color='k',
                   alpha=0.5, edgecolor='none', label=r'$\sigma$')

        ax.annotate(an_txt, [0.05,0.95], xycoords='axes fraction',
                    ha='left', va='top')
        ax.legend(loc='lower right', fancybox=True) #frameon=False)

        ax.set_ylabel(r'Average / Median [$\AA$]')
        ax.set_xlabel(r'X (Along longslit) [pixels]')
        ax.minorticks_on()

        fig.subplots_adjust(left=0.11, right=0.99, bottom=0.12, top=0.96)
        fig.set_size_inches(8,4)
        fig.savefig(pp, format='pdf')

        pp.close()

    if silent == False: mylogger.info('### End residual_wave_cal : '+systime())
#enddef

def cross_check(path, cdir, dbase):
    '''
    Check arc/OH calibration against OH/arc dataset

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
    Created by Chun Ly, 31 May 2018
     - Generate transformed image with iraf.gnirs.nstransform
     - Call mylogger to warn if files exist
     - Change gnirs logfile
     - Move call to iraf.gemini.nsheaders to this function
     - iraf.nstransform does not like suffixes. Using underscores
     - Call OH_check
     - Check for database file before running nsfitcoords and nstransform
    Modified by Chun Ly, 3 June 2018
     - Bug fix: Call OH_check with cross_check=True
     - Call arc_check2 for database_OH case
    Modified by Chun Ly, 8 June 2018
     - Exit when db_file not available
    Modified by Chun Ly, 19 June 2018
     - Call get_database_model
     - Pass function and order to nsfitcoords
    Modified by Chun Ly, 20 June 2018
     - Set nsfitcoords xorder fitting
    '''

    logfile  = path+'QA_wave_cal.log'
    mylogger = glog.log0(logfile)._get_logger()

    timestamp = systime().replace(':','.')
    logfile   = path+'gnirs_'+timestamp+'.log'
    iraf.gemini.gnirs.logfile = logfile

    iraf.gemini.nsheaders("gnirs")

    mylogger.info("Raw data is located in : %s" % path)
    mylogger.info("GNIRS logfile : "+logfile)

    iraf.chdir(path)
    if dbase == 'database/':
        infile  = 'OH_stack.fits'
        outfile = 'fOH_stack_arc.fits'
        lamp    = 'warc_stack.fits'

    if dbase == 'database_OH/':
        infile  = 'arc_stack.fits'
        outfile = 'farc_stack_OH.fits'
        lamp    = 'wOH_stack.fits'

    db_file = dbase+'id'+lamp.replace('.fits','_SCI_1_')

    if not exists(db_file):
        mylogger.warn("Wavelength calibration file not found : "+db_file)
        mylogger.warn("Exiting!!!")
        return
    else:
        if not exists(outfile):
            source0 = 'OH' if '_OH' in dbase else 'arc'
            func0, order0 = get_database_model(path, source0)
            iraf.gnirs.nsfitcoords(infile, outprefix='', outspectra=outfile,
                                   lamp=lamp, database=dbase,
                                   function=func0, lyorder=order0, lxorder=xorder)
        else:
            mylogger.warn('File exists! : '+outfile)

        t_outfile = 't'+outfile
        if not exists(t_outfile):
            iraf.gnirs.nstransform(outfile, outprefix='', outspectra=t_outfile,
                                   database=dbase)
        else:
            mylogger.warn('File exists! : '+t_outfile)

    iraf.chdir(cdir)

    if dbase == 'database/':
        OH_check(path, cross_check=True)

    # + on 03/06/2018
    if dbase == 'database_OH/':
        arc_check2(path, cross_check=True)
#enddef
