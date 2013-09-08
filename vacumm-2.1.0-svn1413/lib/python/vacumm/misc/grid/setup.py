"""Installation script"""

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('grid', parent_package, top_path)
    config.add_extension('_interp_',  sources=['interp.f90'])#,  f2py_options=["--f77exec=gfortran",  "--f90exec=gfortran"])
    config.add_data_files('Makefile')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)

