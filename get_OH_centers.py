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
from IQ_plot import gauss1d

from scipy.optimize import curve_fit

co_dirname = os.path.dirname(__file__)
OH_file = co_dirname+'/rousselot2000.dat'

def group(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last: # Part of the group, bump the end
            last = n
        else: # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last # Yield the last group

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

    i_lines = np.where(OH_spec_mod >= np.max(OH_spec_mod)*0.01)[0]

    lines_set = list(group(i_lines))

    OH_spec_mod_resid = OH_spec_mod.copy()

    rev_lines = np.zeros(len(lines_set))

    for ii in range(len(lines_set)):
        x_avg   = np.int(np.average(lines_set[ii]))
        lam_cen = x0[x_avg]
        sig     = lam_cen / R_spec / (2*np.sqrt(2*np.log(2)))
        p0      = [0.0, OH_spec_mod[x_avg], lam_cen, sig]
        zoom    = np.arange(lines_set[ii][0],lines_set[ii][1])
        bounds = ((-0.001, 0.0, lam_cen-0.5, 0.1),
                  (0.001, 1.25*p0[1], lam_cen+0.5, 1.5*p0[3]))
        popt, pcov = curve_fit(gauss1d, x0[zoom], OH_spec_mod[zoom], p0=p0,
                               bounds=bounds)

        t_mod = gauss1d(x0, *popt)
        OH_spec_mod_resid -= t_mod

        rev_lines[ii] = popt[2]

    if silent == False: log.info('### End main : '+systime())
#enddef

