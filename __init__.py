"""
GNIRSLongSlit

Python 2.7 codes to reduce longslit data from Gemini-N/GNIRS

This code has been tested with X- and J-band spectra with the 111/mm grating

Requirements:
 astropy    >=v1.3
 chun_codes
 matplotlib >= v1.5.3
 numpy      >= v1.11.1

How to use:
import GNIRSLongSlit
"""

import hdr_info
import create_list
import QA_plot

import cleanir
import cleanir_script

reload(hdr_info)
reload(create_list)
reload(QA_plot)

reload(cleanir)
reload(cleanir_script)
