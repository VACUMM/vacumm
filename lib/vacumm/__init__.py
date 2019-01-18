#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2010-2018)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#

from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from warnings import warn
import six
import logging

__project__ = 'vacumm'
__version__ = '3.5.1'
__release__ = '0'
__date__ = '2017-12-12'
__author__ = u'Stéphane Raynaud, Jonathan Wilkins, Guillaume Charria'
__email__ = 'raynaud@actimar.fr, wilkins@actimar.fr, charria@ifremer.fr'
__copyright__ = 'Copyright (c) 2010-2016 Actimar/IFREMER'
__description__ = """
VACUMM library

A collection of tools for VACUMM python codes.

Content:

    - misc: generic library
    - data: data management
    - diag: advanced diagnostics
    - tide: tidal tools
    - sphinxext: extensions to sphinx doc generator
    - bathy: bathymetric tools
    - data: library to read observations and modelled fields
    - validator: tools for model validation
    - sphinxext: sphinx extensions
    - report: reporting
    - markup: tools for html page generation

"""

docfiller_verbose = False


def help(text=None, recent=False):
    """Open VACUMM website in a web browser and optionally search for a string

    :Params:

        - **text**, optional: String to search for.
        - **recent**, optional: Use the most recent version of the documentation.
    """
    key = 'url_recent' if recent else 'url'
    url = VACUMM_CFG['vacumm'][key]
    if url is None:
        url = 'http://www.ifremer.fr/vacumm'
    from webbrowser import open
    if text is not None:
        if not isinstance(text, six.string_types):
            if hasattr(text, 'func_name'):
                text = text.__name__
            else:
                text = text.__class__.__name__
        if not text.startswith('/'):
            text = '/search.html?q=%s&check_keywords=yes&area=default' % text
        url += text
    open(url, new=2)


class VACUMMError(Exception):
    """Standard VACUMM error (exception)"""


class VACUMMWarning(UserWarning):
    """Standard VACUMM warning"""


def vacumm_warn(message, stacklevel=2):
    """Issue a :class:`VACUMMWarning`"""
    warn(message, VACUMMWarning, stacklevel=stacklevel)


vcwarn = vacumm_warning = vacumm_warn


#: Current VACUMM configuration (default+common+user)
VACUMM_CFG = CFG = None
from . import config
config.load_cfg('com')
config.load_cfg('user')


# Logging
class NullHandler(logging.Handler):
    def emit(self, record):
        pass
logging.getLogger("vacumm").addHandler(NullHandler())
