# -*- coding: utf8 -*-
"""Installation script"""
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    import os
    config = Configuration('satellite', parent_package, top_path)
    return config

