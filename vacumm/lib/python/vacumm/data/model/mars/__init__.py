"""
Data_reader tools

Data_reader submodule can be used with vacumm.data_reader.*
"""

import get_ftp
import get_opendap
import get_cp
import os as _os, locale as _locale
_os.environ['LC_NUMERIC'] = 'en_US.UTF-8'
_locale.setlocale(_locale.LC_NUMERIC, 'en_US.UTF-8')

#from mars3D import MARS3D
