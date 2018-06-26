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

from fit import gaussian, gaussian_R
from IQ_plot import gauss1d

from scipy.optimize import curve_fit

import glog

co_dirname = os.path.dirname(__file__)
OH_file = co_dirname+'/rousselot2000.dat'

bbox_props = dict(boxstyle="square,pad=0.15", fc="w", alpha=0.75, ec="none")

n_multi = 8

def gauss_multi(x, a0, a1, a2, a3, a4, a5, a6, a7,
                l0, l1, l2, l3, l4, l5, l6, l7,
                s0, s1, s2, s3, s4, s5, s6, s7):
    m0 = a0 * np.exp(-(x - l0)**2 / (2 * s0**2))
    m1 = a1 * np.exp(-(x - l1)**2 / (2 * s1**2))
    m2 = a2 * np.exp(-(x - l2)**2 / (2 * s2**2))
    m3 = a3 * np.exp(-(x - l3)**2 / (2 * s3**2))
    m4 = a4 * np.exp(-(x - l4)**2 / (2 * s4**2))
    m5 = a5 * np.exp(-(x - l5)**2 / (2 * s5**2))
    m6 = a6 * np.exp(-(x - l6)**2 / (2 * s6**2))
    m7 = a7 * np.exp(-(x - l7)**2 / (2 * s7**2))

    return m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7

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

    Modified by Chun Ly, 13 June 2018
     - Fix rev_lines and rev_int (wrong indexing)
     - Switch to micron units; Do subplots_adjust for white space
     - Plot aesthetics: label plots, change page size, legend
     - Include parameters/assumptions in title
     - Plot aesthetics: Vertical lines for all possible OH skylines
     - Plot aesthetics: Label OH skylines
     - Plot aesthetics: Limit vertical lines for OH skylines
     - Plot aesthetics: group and label OH skylines to avoid overlap
     - Group lines (<5 Ang) before annotation
     - Opaque white background behind lines
     - Plot aesthetics: Draw OH lines to top of subplots
     - Plot aesthetics: Remove legend; adjust white space
     - WARN if more than six lines
     - Attempt to constraint fit using bounds but that did not work
     - Fix case if curve_fit solutions are outside of spectral range
    Modified by Chun Ly, 21 June 2018
     - Write npz file containing use lines that are grouped together
     - mylogger call for writing files
     - Switching back to uncompressed npz - Got interpret file pickle error (IOError)
    Modified by Chun Ly, 25 June 2018
     - Bug fix: Crash on call to group_OH_lines. Require in_zoom to not be empty array
    '''
    
    # + on 09/01/2018
    logfile  = rawdir+'get_OH_centers.log'
    mylogger = glog.log0(logfile)._get_logger()

    if silent == False: mylogger.info('### Begin main : '+systime())

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

    y_max = max(OH_spec_mod) * 1.25

    i_lines = np.where(OH_spec_mod >= np.max(OH_spec_mod)*0.01)[0]

    lines_set = list(group(i_lines))

    OH_spec_mod_resid = OH_spec_mod.copy()

    rev_lines = [] #np.zeros(len(lines_set))
    rev_int   = [] #np.zeros(len(lines_set))

    nrows = 3
    fig, ax = plt.subplots(nrows=nrows)

    use_lines = [] # + on 21/06/2018

    xlim_arr = []
    dx = (x0[-1]-x0[0])/3.0
    for aa in range(nrows):
        ax[aa].plot(x0/1e4, OH_spec_mod, color='black', zorder=3,
                    label="Rousselot (2000)")
        xlim = np.array([x_min+dx*aa,x_min+dx*(aa+1)])
        xlim_arr.append(xlim)

        ax[aa].set_xlim(xlim/1e4)
        ax[aa].set_ylim([-10,y_max])

        # Draw vertical lines for all possible OH skylines
        for val in OH_lines[in_rge]:
            ax[aa].axvline(x=val/1e4, color='black', linewidth=0.25,
                           linestyle=':', zorder=2)

    for ii in range(len(lines_set)):
        x_avg   = np.int(np.average(lines_set[ii]))
        tl_min, tl_max = x0[lines_set[ii][0]], x0[lines_set[ii][1]]
        in_zoom = np.where((OH_lines >= tl_min) & (OH_lines <= tl_max))[0]

        if len(in_zoom) > 0:
            group_lines, group_int = group_OH_lines(OH_lines[in_zoom],
                                                    OH_int[in_zoom])

            # print ii, group_lines

            sig = group_lines / R_spec / (2*np.sqrt(2*np.log(2)))
            peak0 = OH_spec_mod[np.int_((group_lines-x_min)/(x0[1]-x0[0]))]

            if len(group_lines) == 1:
                p0 = [0.0, peak0, group_lines[0], sig[0]]
                #plt.axvline(group_lines[0]/1e4, color='blue')
            else:
                if len(group_lines) > n_multi:
                    log.warn('More than 8 lines found, N='+str(len(group_lines))+'!!!')
                t_peak0 = peak0.tolist()
                t_lines = group_lines.tolist()
                t_sig   = sig.tolist()

                t_peak0 += np.zeros(n_multi-len(group_lines)).tolist()
                t_lines += np.zeros(n_multi-len(group_lines)).tolist()
                t_sig   += np.zeros(n_multi-len(group_lines)).tolist()

                p0 = t_peak0
                p0 += t_lines
                p0 += t_sig

                p0 = np.array(p0)

                zoom    = np.arange(lines_set[ii][0],lines_set[ii][1])
                # print p0
                #bounds = ((-0.001, 0.0, lam_cen-0.5, 0.1),
                #          (0.001, 1.25*p0[1], lam_cen+0.5, 1.5*p0[3]))

            if len(group_lines) == 1:
                popt, pcov = curve_fit(gauss1d, x0[zoom], OH_spec_mod[zoom],
                                       p0=p0)
                t_mod = gauss1d(x0, *popt)

                use_lines.append([popt[2]]) # + on 21/06/2018

                rev_lines.append(popt[2])
                rev_int.append(popt[1])

                i_ax = [xx for xx in range(nrows) if
                        (popt[2] >= xlim_arr[xx][0] and popt[2] <= xlim_arr[xx][1])][0]
                ax[i_ax].annotate('%.2f' % popt[2], [popt[2]/1e4, y_max*0.99], ha='center',
                                  va='top', rotation=90, fontsize=4, bbox=bbox_props)
            else:
                #low_bound = tuple([0] * n_multi) + tuple(p0[n_multi:2*n_multi]-0.5) + \
                    #            tuple([0] * n_multi)
                #up_bound  = tuple(p0[0:n_multi]*1.25+0.1) + tuple(p0[n_multi:2*n_multi]+0.5) + \
                    #            tuple(p0[2*n_multi:]+1)
                popt, pcov = curve_fit(gauss_multi, x0[zoom], OH_spec_mod[zoom], p0=p0)
                t_mod = gauss_multi(x0, *popt)

                t_loc = popt[n_multi:2*n_multi] # Line wavelengths (Ang)
                t_str = popt[0:n_multi]         # Line peak strength

                nonzero = np.where(t_loc != 0)[0]
                use_lines.append(t_loc[nonzero].tolist()) # + on 21/06/2018

                wave0 = t_loc[nonzero]
                wave0.sort()

                # Check that lines are within range
                in_rge = np.where((wave0 >= x_min) & (wave0 <= x_max))[0]
                wave0 = wave0[in_rge]
                t_str = t_str[nonzero[in_rge]]

                rev_lines += wave0.tolist()
                rev_int   += t_str.tolist()

                wave0 = np.array(wave0)
                # Label lines
                skip = np.zeros(len(wave0))
                for ww in range(len(wave0)):
                    if skip[ww] == 0:
                        w_diff = wave0[ww:]-wave0[ww]
                        t_close = np.where(w_diff <= 5)[0]
                        close = np.arange(ww,len(wave0))[t_close]
                        str_comb = "\n".join(['%.2f' % val for
                                              val in wave0[close]])
                        w_cen = np.average(wave0[close])
                        #0.5*(wave0[0]+wave0[-1]) #np.average(wave0)

                        i_ax = [xx for xx in range(nrows) if
                                (w_cen >= xlim_arr[xx][0] and w_cen <= xlim_arr[xx][1])][0]
                        ax[i_ax].annotate(str_comb, [w_cen/1e4, y_max*0.99], ha='center',
                                          va='top', rotation=90, fontsize=4, bbox=bbox_props)
                        skip[close] = 1
                    #endif
                #endfor

        # print '## t_mod : ', np.min(t_mod), np.max(t_mod)

        OH_spec_mod_resid -= t_mod

    for aa in range(nrows):
        ax[aa].plot(x0/1e4, OH_spec_mod_resid, linestyle='dashed', color='blue',
                    zorder=3, label='Residual')

        for t_line in rev_lines:
            ax[aa].axvline(x=t_line/1e4, color='red', linestyle='--',
                           linewidth=1.0, zorder=1)

    l_tab = Table([rev_lines, rev_int])
    out_file = rawdir+'rousselot2000_convl.dat'
    mylogger.info('Writing : '+out_file)
    asc.write(l_tab, out_file, format='no_header')

    ann_txt  = r'%.2f $\mu$m; %s; ' % (gratwave, o_tab0['filter2'])
    ann_txt += '%s; %.3f" slit; R = %i' % (o_tab0['grating'], slitwidth, R_spec)
    ax[0].set_title(ann_txt)
    #leg = ax[0].legend(loc='upper right', fancybox=True) #, frameon=False)
    #leg.get_frame().set_alpha(0.75)

    ax[2].set_xlabel(r'Wavelength [$\mu$m]')
    plt.subplots_adjust(left=0.08, right=0.95, bottom=0.06, top=0.96,
                        hspace=0.12)

    fig.set_size_inches(6,8)
    out_pdf = out_file.replace('.dat','.pdf')
    mylogger.info('Writing : '+out_pdf)
    fig.savefig(out_pdf)

    # Write npz file containing final grouping | + on 21/06/2018
    savez_file = out_file.replace('.dat','.npz')
    mylogger.info('Writing : '+savez_file)
    np.savez(savez_file, use_lines=use_lines)

    if silent == False: mylogger.info('### End main : '+systime())
#enddef

