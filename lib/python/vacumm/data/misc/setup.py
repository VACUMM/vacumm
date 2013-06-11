"""Installation script"""
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    import os
    config = Configuration('misc', parent_package, top_path)
    config.add_data_files('coloc.ini')
    config.add_data_files('dataset.ini')
    config.add_data_files('profile.ini')
    config.add_data_files('satellite.ini')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)

