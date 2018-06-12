"""
get_OH_centers
==============

Convolves Rousselot model with lower resolution to determine OH night skylines
center
"""

import sys, os

from chun_codes import systime

from os.path import exists
import commands
from astropy.io import ascii as asc

import numpy as np
from math import pi

import matplotlib.pyplot as plt

from astropy.table import Table
from astropy import log

import glob

from OH_stack import gaussian, gaussian_R

co_dirname = os.path.dirname(__file__)
OH_file = co_dirname+'/rousselot2000.dat'

def main(rawdir, silent=False, verbose=True):

    '''
    Main function for get_OH_centers

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
    Created by Chun Ly, 12 June 2018
    '''
    
    if silent == False: log.info('### Begin main : '+systime())

    if exists(OH_file):
        if silent == False: mylogger.info('Reading : '+OH_file)
        OH_data  = np.loadtxt(OH_file)
        OH_lines = OH_data[:,0]
        OH_int   = OH_data[:,1]

    infile = rawdir + 'hdr_info.QA.tbl'
    tab0  = asc.read(infile, format='fixed_width_two_line')
    i_obj = [xx for xx in range(len(tab0)) if tab0['QA'][xx] == 'obj'][0]
    o_tab0 = tab0[i_obj]
    gratwave = o_tab0['gratwave']

    slitwidth = np.float(o_tab0['slit'].split('arcsec')[0])
    if '111/mm' in o_tab0['grating']:
        if o_tab0['filter2'] == 'X_G0518':
            R_spec = 6600.0/(slitwidth/0.3)
            x_diff = 0.06
        if o_tab0['filter2'] == 'J_G0517':
            R_spec = 7200.0/(slitwidth/0.3)
            x_diff = 0.07
    
    x_min = (gratwave-x_diff) * 1.E4
    x_max = (gratwave+x_diff) * 1.E4
    npix  = np.round((x_max-x_min)/0.15)

    x0 = x_min + 0.15*np.arange(npix)

    OH_spec_mod = np.zeros(len(x0))
    in_rge = np.where((OH_lines >= x0[0]) & (OH_lines <= x0[-1]))[0]
    for ii in range(len(in_rge)):
        idx = in_rge[ii]
        temp = OH_int[idx] * gaussian_R(x0, OH_lines[idx], R_spec)
        OH_spec_mod += temp


    if silent == False: log.info('### End main : '+systime())
#enddef

