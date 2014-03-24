"""Installation script"""
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('misc', parent_package, top_path)
#    config.add_data_files('logo_vacumm.gif')
    config.add_subpackage('grid')
    config.add_subpackage('phys')
    config.add_subpackage('axml')
#   config.add_subpackage('easyPypar')
    config.add_data_files('config.cfg')
    config.add_data_files('cpt/*')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)

