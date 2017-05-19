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

from astropy.io import fits
from astropy import log

def arc_check(path, arcs='', silent=False, verbose=True):

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

    for nn in xrange(n_arcs):
        print rnc_files[nn]
        hdu0 = fits.open(rnc_files[nn])
        im0  = hdu0['SCI'].data
        f    = open(d_files[nn], 'r')
        str0 = f.readlines()

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
                l_val0[cc,rr] = np.float(temp1[2])
                l_val[cc,rr]  = np.float(temp1[3])

        #print len(im0[0])
        plt.imshow(im0, cmap='Greys', origin='lower')
    #endfor

    if silent == False: log.info('### End QA_wave_cal : '+systime())
#enddef

