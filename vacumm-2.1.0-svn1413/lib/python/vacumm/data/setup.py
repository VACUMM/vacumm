"""Installation script"""
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    import os
    config = Configuration('data', parent_package, top_path)
    config.add_subpackage('model')
    config.add_subpackage('satellite')
    config.add_subpackage('misc')
    config.add_subpackage('in_situ')
    return config

