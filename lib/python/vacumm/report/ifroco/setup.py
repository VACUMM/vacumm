"""Installation script"""
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('ifroco', parent_package, top_path)
    config.add_data_files('config.cfg')
    return config
