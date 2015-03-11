# -*- coding: utf8 -*-
"""Installation script"""
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    import os
    config = Configuration('model', parent_package, top_path)
    config.add_subpackage('mars')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)

