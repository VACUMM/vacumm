"""Installation script"""
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('vacumm', parent_package, top_path)
    config.add_subpackage('data')
    config.add_subpackage('misc')
    config.add_subpackage('diag')
    config.add_subpackage('sphinxext')
    config.add_subpackage('markup')
    config.add_subpackage('bathy')
    config.add_subpackage('tide')
    config.add_subpackage('validator')
    config.add_data_files('config.cfg')
    return config
