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
from astropy import log

def arc_check(path, arcs='', out_pdf='', silent=False, verbose=True):

    '''
    Generate plot illustrating wavelength calibration from IRAF database

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
    Created by Chun Ly, 18 May 2017
    Modified by Chun Ly, 19 May 2017
     - Fix bug with scatter plot (wrong indexing)
     - Fix bug with l_val0 and l_val (those were flipped)
     - Draw dashed lines connecting the points
     - Adjust margins, add labels
    '''
    
    if silent == False: log.info('### Begin QA_wave_cal : '+systime())

    if arcs == '':
        arc_list = path + 'arc.lis'
        if silent == False: log.info('### Reading : '+arc_list)
        arcs = np.loadtxt(arc_list, dtype=type(str))

    n_arcs = len(arcs)

    rnc_files = [path+'rnc'+file0 for file0 in arcs]
    d_files   = [path+'database/idwrnc'+file0.replace('.fits','_SCI_1_') \
                 for file0 in arcs]

    out_pdf = path+'arc_check.pdf' if out_pdf == '' else path+out_pdf
    pp = PdfPages(out_pdf)

    for nn in xrange(n_arcs):
        hdu0 = fits.open(rnc_files[nn])
        im0  = hdu0['SCI'].data

        fig, ax = plt.subplots()
        ax.imshow(im0, cmap='Greys', origin='lower')

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
            temp = temp.replace('wrnc'+arcs[nn].replace('.fits','[SCI,1]['),'')
            x_cols[cc] = np.int(temp.split(',')[0])

            for rr in xrange(n_features[cc]): #beg_mark[cc]+1, func_mark[cc]):
                temp1 = str0[feat_mark[cc]+1+rr].replace('       ','')
                temp1 = temp1.split(' ')
                temp1 = [val for val in temp1 if val != '']
                y_val[cc,rr]  = np.float(temp1[1])
                l_val[cc,rr]  = np.float(temp1[2])
                l_val0[cc,rr] = np.float(temp1[3])
            #endfor
            x_temp = np.repeat(x_cols[cc],n_features[cc])
            y_temp = y_val[cc,0:n_features[cc]]
            ax.scatter(x_temp, y_temp, facecolor='none', edgecolor='red',
                       marker='o', s=25)
        #endfor

        line_list = list(set(l_val0.reshape(n_cols*n_max).tolist()))
        print line_list
        for ll in range(len(line_list)):
            if line_list[ll] != 0:
                #print line_list[ll]
                x_idx, y_idx = np.where(l_val0 == line_list[ll])
                sort0 = np.argsort(x_cols[x_idx])
                ax.plot(x_cols[x_idx][sort0], y_val[x_idx,y_idx][sort0], 'r--')
        ax.set_xlabel('X [pixels]')
        ax.set_ylabel('Y [pixels]')
        ax.minorticks_on()
        plt.subplots_adjust(left=0.10, bottom=0.025, top=0.99,
                           right=0.99, wspace=0.02, hspace=0.02)
        fig.set_size_inches(8,11)
        fig.savefig(pp, format='pdf')
    #endfor

    pp.close()
    if silent == False: log.info('### End QA_wave_cal : '+systime())
#enddef
