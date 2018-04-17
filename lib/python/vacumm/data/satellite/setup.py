# -*- coding: utf8 -*-
"""Installation script"""
from __future__ import absolute_import
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    import os
    config = Configuration('satellite', parent_package, top_path)
    return config

