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

def gauss_six(x, a0, a1, a2, a3, a4, a5, l0, l1, l2, l3, l4, l5,
              s0, s1, s2, s3, s4, s5):
    m0 = a0 * np.exp(-(x - l0)**2 / (2 * s0**2))
    m1 = a1 * np.exp(-(x - l1)**2 / (2 * s1**2))
    m2 = a2 * np.exp(-(x - l2)**2 / (2 * s2**2))
    m3 = a3 * np.exp(-(x - l3)**2 / (2 * s3**2))
    m4 = a4 * np.exp(-(x - l4)**2 / (2 * s4**2))
    m5 = a5 * np.exp(-(x - l5)**2 / (2 * s5**2))

    return m0 + m1 + m2 + m3 + m4 + m5

def group(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last: # Part of the group, bump the end
            last = n
        else: # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last # Yield the last group

def group_OH_lines(OH_lines, OH_int):
    str_uniq = np.array(list(set(OH_int)))
    str_uniq.sort()

    bright = np.where(str_uniq >= 0.05*max(str_uniq))[0]

    best_lines = np.zeros(len(bright))
    for ii in range(len(bright)):
        idx = np.where(OH_int == str_uniq[bright[ii]])[0]
        best_lines[ii] = np.sum(OH_lines[idx]*OH_int[idx])/np.sum(OH_int[idx])
    #endfor
    return best_lines, str_uniq[bright]
#enddef


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

    rev_lines = [] #np.zeros(len(lines_set))
    rev_int   = [] #np.zeros(len(lines_set))

    for ii in range(len(lines_set)):
        x_avg   = np.int(np.average(lines_set[ii]))
        tl_min, tl_max = x0[lines_set[ii][0]], x0[lines_set[ii][1]]
        in_zoom = np.where((OH_lines >= tl_min) & (OH_lines <= tl_max))[0]

        group_lines, group_int = group_OH_lines(OH_lines[in_zoom],
                                                OH_int[in_zoom])

        # print ii, group_lines

        sig = group_lines / R_spec / (2*np.sqrt(2*np.log(2)))
        peak0 = OH_spec_mod[np.int_((group_lines-x_min)/(x0[1]-x0[0]))]

        if len(group_lines) == 1:
            p0 = [0.0, peak0, group_lines[0], sig[0]]
            plt.axvline(group_lines[0], color='blue')
        else:
            t_peak0 = peak0.tolist()
            t_lines = group_lines.tolist()
            t_sig   = sig.tolist()

            t_peak0 += np.zeros(6-len(group_lines)).tolist()
            t_lines += np.zeros(6-len(group_lines)).tolist()
            t_sig   += np.zeros(6-len(group_lines)).tolist()

            p0 = t_peak0
            p0 += t_lines
            p0 += t_sig


        zoom    = np.arange(lines_set[ii][0],lines_set[ii][1])
        # print p0
        #bounds = ((-0.001, 0.0, lam_cen-0.5, 0.1),
        #          (0.001, 1.25*p0[1], lam_cen+0.5, 1.5*p0[3]))
        if len(group_lines) == 1:
            popt, pcov = curve_fit(gauss1d, x0[zoom], OH_spec_mod[zoom],
                                   p0=p0)
            t_mod = gauss1d(x0, *popt)
            rev_lines.append(popt[2])
            rev_int.append(popt[1])
        else:
            popt, pcov = curve_fit(gauss_six, x0[zoom], OH_spec_mod[zoom],
                                   p0=p0)
            t_mod = gauss_six(x0, *popt)

            rev_lines += popt[range(len(group_lines),2*len(group_lines))]
            rev_int   += popt[0:len(group_lines)]

        # print '## t_mod : ', np.min(t_mod), np.max(t_mod)

        OH_spec_mod_resid -= t_mod

        l_tab = Table([rev_lines, rev_int])
        asc.write(l_tab, rawdir+'rousselot2000_convl.dat', format='no_header')
    if silent == False: log.info('### End main : '+systime())
#enddef

