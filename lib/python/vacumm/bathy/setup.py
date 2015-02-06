"""Installation script"""

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('bathy', parent_package, top_path)
    config.add_data_files('vacumm.cfg')
    config.add_data_files('bathy.gridded.cfg')
#    config.add_data_files('bathy.xyz.default.cfg')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)

