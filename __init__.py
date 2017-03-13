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

# + on 12/03/2017
gnirs_2017a = ['DEEP05', 'DEEP06', 'DEEP07', 'Keck03', 'Keck27', 'MMT37']

__all__ = [gnirs_2017a]

import hdr_info
import create_list
import QA_plot

import cleanir
import cleanir_script
import IQ_plot

reload(hdr_info)
reload(create_list)
reload(QA_plot)

reload(cleanir)
reload(cleanir_script)
reload(IQ_plot)
