# -*- coding: utf-8 -*-
# Copyright or © or Copr. Actimar/IFREMER (2010-2018)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
# Inits
import os, shutil
rootdir = os.path.dirname(__file__)
#sys.path.insert(0, os.path.abspath(os.path.join(rootdir, 'doc/sphinx/source')))

# Revision number
# - from __init__.py
f = open(os.path.join(rootdir, 'lib/vacumm/__init__.py'))
for line in f:
    line = line[:-1].strip()
    if line.startswith('__version__'):
        exec(line)
        release = __version__
        break
f.close()
version_sphinx = release
release_sphinx = release

# Infos
name = "vacumm"
version = release
description = 'A library for ocean science'
long_description = ("A library and a collection of scripts "
                    "for ocean science, mainly designed for "
                    "data analysis and model validation")
author = 'Actimar / IFREMER'
author_email = 'raynaud@actimar.fr or charria@ifremer.fr'
maintainer = "Actimar / IFREMER"
maintainer_email = author_email
plateform = 'Linux/UNIX/BSD'
license = "CeCiLL"
url = "http://www.ifremer.fr/vacumm"
classifiers = ["Development Status :: 4 - Beta",
               "Intended Audience :: Science/Research",
               "License :: CeCiLL",
               "Programming Language :: Python :: 2",
               "Topic :: Scientific/Engineering :: GIS",
               "Topic :: Scientific/Engineering :: Physics",
               "Topic :: Scientific/Engineering :: Mathematics",
               "Topic :: Scientific/Engineering :: Atmospheric Science",
               "Topic :: Software Development :: Libraries :: Python Modules",
               "Operating System :: POSIX",
               "Operating System :: UNIX",
               "Operating System :: MacOS :: MacOS X",
               ]

# Setup
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from glob import glob
def configuration(parent_package='',top_path=None):

    # Initialize
    config = Configuration(None,
        parent_package=parent_package,
        top_path=top_path,
    )

    # Set options
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        quiet=True,
        delegate_options_to_subpackages=True,
    )

    # Add bin scripts
    scripts = []
    for pat in ['*.py', '*.bash', '*.sh']:
        pat = os.path.join('bin', pat)
        scripts.extend(glob(pat))
    config.add_scripts(*scripts)

    # Add special files
    config.add_data_files(('vacumm', ['LICENSE', 'CHANGES.rst']))
    config.add_subpackage('vacumm', 'lib/vacumm')
    config.add_data_dir(('vacumm/vacumm-config', 'config')) # not needed
    config.add_data_dir(('vacumm/vacumm-data', 'data')) # all data files
    config.add_data_dir(('vacumm/vacumm-scripts/test', 'scripts/test')) # test scripts
    config.add_data_dir(('vacumm/vacumm-scripts/tutorials', 'scripts/tutorials')) # tutorials
    config.add_data_dir(('vacumm/vacumm-scripts/courses', 'scripts/courses')) # courses

    # Add extensions
    config.add_extension('vacumm.fortran.interp', ['src/interp.f90'])

    return config

# Add the --cfgfile option to install command
from numpy.distutils.command.install_data import install_data as numpy_install_data
#from numpy.distutils.command.sdist import sdist as numpy_sdist

class vacumm_install_data(numpy_install_data):
    user_options = numpy_install_data.user_options + [
        ('cfgfiles=', 'c', "VACUMM main and secondary "
            "configuration files to replace defaut values")]
    def initialize_options (self):
        numpy_install_data.initialize_options(self)
        self.cfgfiles = (hasattr(self.distribution, 'cfgfiles')
            and self.distribution.cfgfiles or None)
    def finalize_options(self):
        numpy_install_data.finalize_options(self)
        self.distribution.cfgfiles = self.cfgfiles
        if self.cfgfiles is not None:
            cfgfiles = []
            for path in self.cfgfiles.split(','):
#                self.add_data_files('vacumm/vacumm-config', path)
                path = os.path.expanduser(path)
                cfgfiles.extend(glob(os.path.abspath(path)))
            self.data_files.append(('vacumm/vacumm-config', cfgfiles)) # FIXME: CFG INSTALL DIR NOT CLEAN

from numpy.distutils.command.install import install as numpy_install
class vacumm_install(numpy_install):
    user_options = numpy_install.user_options + [
        ('cfgfiles=', None, "VACUMM main and secondary "
            "configuration files to replace defaut values"),
    ]
    def initialize_options (self):
        numpy_install.initialize_options(self)
        self.cfgfiles = None
    def finalize_options(self):
        numpy_install.finalize_options(self)
        self.distribution.cfgfiles = self.cfgfiles

import distutils.command.bdist_rpm
class vacumm_bdist_rpm(distutils.command.bdist_rpm.bdist_rpm):
    def initialize_options (self, *args, **kwargs):
        distutils.command.bdist_rpm.bdist_rpm.initialize_options(self)
        import vacumm
        self.distribution_name = name
        self.group = 'Utilities'
        self.release = vacumm.__release__
        self.vendor = vacumm.__copyright__
        self.requires = [
            'python >= 2.6',
            'python-configobj >= 4.6',
            'cdat-lite >= 6.0rc2',
            'sphinx-fortran >= 1.0',
        ]
        #self.install_script = 'rpm-install'
        #self.post_install = 'rpm-post-install'
        #self.pre_uninstall = 'rpm-pre-uninstall'

if __name__ == '__main__':

    # Setup config file
    if not os.path.exists('setup.cfg'):
        shutil.copy('setup.cfg.simple', 'setup.cfg')

    # Lauch setup
    setup(name='vacumm',
          version=version,
          description=description,
          long_description=long_description,
          author=author,
          author_email=author_email,
          maintainer=author,
          maintainer_email=author_email,
          license=license,
          url=url,
          package_dir= {'':'lib'},
          py_modules=['vcmq', 'vcmqm'],
          classifiers=classifiers,
          cmdclass={'install':vacumm_install,
                    'install_data':vacumm_install_data,
                    'bdist_rpm':vacumm_bdist_rpm},
          configuration=configuration,
#        ext_modules=[ext_interp]

    )



