"""
GNIRSLongSlit

Python 2.7 codes to reduce longslit data from Gemini-N/GNIRS

This code has been tested with X- and J-band spectra with the 111/mm grating

Requirements:
 aplpy      >= v1.1.1
 astropy    >= v1.3
 ccdproc    >= v1.2.0
 chun_codes (https://github.com/astrochun/chun_codes)
 matplotlib >= v1.5.3
 scipy      >= 0.18.1
 numpy      >= v1.11.1
 pyraf      >= 2.2.dev
  - Install through astroconda to get gemini IRAF package included

How to use:
import GNIRSLongSlit as gnirsls
 - This will import most modules.

For the main reduction module, reduce, it needs to load PyRAF which delays
the import. Thus, a separate import is needed:
from GNIRSLongSlit import reduce
"""

# + on 12/03/2017
gnirs_2017a = ['DEEP05', 'DEEP06', 'DEEP07', 'DEEP10', 'DEEP14', 'DEEP15',
               'DEEP17', 'DEEP18', 'DEEP19', 'Keck03', 'Keck14', 'Keck27',
               'MMT08', 'MMT37']

gnirs_2017b = ['DEEP26', 'DEEP28']

__all__ = [gnirs_2017a, gnirs_2017b]

import hdr_info
import create_list
import create_obs_summary_table

import cleanir_script

import QA_plot
import IQ_plot
import align_check

import symlink
import dir_check

import iraf_get_subset
import file_handling
#import QA_wave_cal
import remove_bias_level

reload(hdr_info)
reload(create_list)
reload(create_obs_summary_table)

reload(cleanir_script)

reload(QA_plot)
reload(IQ_plot)
reload(align_check)

reload(symlink)
reload(dir_check)

reload(iraf_get_subset)
reload(file_handling)
#reload(QA_wave_cal)
reload(remove_bias_level)
