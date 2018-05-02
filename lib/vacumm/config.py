# -*- coding: utf8 -*-
"""Utilities to manage VACUMM basic configuration of modules"""
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


from __future__ import absolute_import
from __future__ import print_function
import os, sys
import subprocess,  shlex
from warnings import warn
import six.moves.urllib.request, six.moves.urllib.parse, six.moves.urllib.error
import six
from six.moves import input

from configobj import ConfigObj
from validate import Validator, VdtMissingValue

import vacumm
from vacumm import VACUMMError

__all__ = ['VACUMM_CFGSPECS', 'VACUMM_CFGSPECS_FILE', 'VACUMM_VDT', 'VACUMM_VDT_LIVE',
    'VCValidator', 'check_data_file', 'download_file', 'edit_file', 'get_cfg_checked',
    'get_cfg_path', 'get_com_conf_dir', 'get_com_conf_file', 'get_data_dir',
    'get_default_config', 'get_dir_dict', 'get_dist_dir', 'get_dl_dir',
    'get_home_conf_dir', 'get_scripts_dir', 'get_tut_dir', 'get_user_conf_file',
    'get_user_conf_dir', 'load_cfg', 'save_config_value', 'set_config_value',
    'edit_user_conf_file', 'print_config', 'get_config_value', 'data_sample',
    'write_default_config', 'get_lib_dir', 
    ]

def get_lib_dir():
    """Directory of the :mod:`vacumm` package"""
    return os.path.abspath(os.path.dirname(__file__))


def get_dist_dir():
    """Upper directory of VACUMM distribution tree or ``None`` for an installed package"""
    lib_dir = get_lib_dir()
    dist_dir = os.path.abspath(os.path.join(lib_dir, '../..'))
    if not os.path.exists(os.path.join(dist_dir, 'setup.py')): return
    return dist_dir

def get_data_dir(raiseerr=False):
    """Get the data directory absolute path

    This directory contains data samples and other needed files.
    It can be at two different places, depending on if the library is
    an installed version or a developers version.

    - If installed : in :file:`vacumm-data` subdirectory in the
      installed package directory (see :meth:`get_lib_dir`).
    - Else in the :file:`data` subdirectory of the main distribution tree
      (see :meth:`get_dist_dir`).

    .. warning:: It raises an :class:`VACUMMError` error if not found
        and if ``raiseerr`` is True, else return ''.

    """
    # Installed librairy
    lib_dir = get_lib_dir()
    data_dir = os.path.join(lib_dir, 'vacumm-data')
    if os.path.exists(data_dir):
        return data_dir

    # Distributed library (dev)
    dist_dir = get_dist_dir()
    if dist_dir is not None:
        data_dir = os.path.join(dist_dir, 'data')
        if os.path.exists(data_dir):
            return data_dir

    if raiseerr: raise VACUMMError("Can't find a valid data directory")
    return ''

def data_sample(data_file, mode='all'):
    """Transform the relative path of a sample file to an absolute path

    This can be either in the VACUMM data directory (:func:`get_data_dir`),
    or in the UVCDAT data directory.

    :Params:

        - **data_file**: File basename.
        - **mode**, optional:

            - ``"all"`` or ``None``: Look first into VACUMM dir then
              and into UVCDAT dir if not found. If finally not found,
              return the VACUMM sample path.
            - ``"vacumm"``: Look only into the VACUMM data dir only.
            - ``"uvcdat"``: Look only into the UVCDAT data dir only.

    """
    if mode is None: mode = 'all'
    mode = str(mode)

    # VACUMM
    if mode.startswith('a') or not mode.startswith('u'):
        path = vacpath = os.path.join(get_data_dir(raiseerr=True), data_file)
        if os.path.exists(path) or not mode.startswith('a'): return path

    # UVCDAT
    path = os.path.join(sys.prefix, 'sample_data', data_file)
    if os.path.exists(path) or mode.startswith('u'): return path
    return vacpath


def get_scripts_dir(subdir=None, raiseerr=False):
    """Get the scripts directory absolute path

    This directory contains examples of script.
    It can be at two different places, depending on if the library is
    an installed version or a developers version.

    - If installed : in :file:`vacumm-scripts` subdirectory in the
      installed package directory (see :meth:`get_lib_dir`).
    - Else in the :file:`scripts` subdirectory of the main distribution tree
      (see :meth:`get_dist_dir`).

    .. warning:: It raises an :class:`VACUMMError` error if not found
        and if ``raiseerr`` is True.
    """
    # Installed librairy
    lib_dir = get_lib_dir()
    scripts_dir = os.path.join(lib_dir, 'vacumm-scripts')
    if not os.path.exists(scripts_dir): scripts_dir = None

    # Distributed library (dev)
    if scripts_dir is None:
        dist_dir = get_dist_dir()
        if dist_dir is not None:
            scripts_dir = os.path.join(dist_dir, 'scripts')
            if not os.path.exists(scripts_dir): scripts_dir = None

    # Not found
    if scripts_dir is None:
        if raiseerr: raise VACUMMError("Can't find a valid scripts directory")
        return ''

    # Subdir
    if subdir and isinstance(subdir, six.string_types):
        scripts_dir = os.path.join(scripts_dir, subdir)
        if not os.path.exists(scripts_dir):
            raise VACUMMError("Invalid subdirectory of the scripts directory: "+subdir)

    return scripts_dir


def get_tut_dir(raiseerr=False):
    """Get the tutorials directory absolute path

    This directory contains all tutorials used by the documentation.
    It can be at two different places, depending on if the library is
    an installed version or developers version.

    - If installed : in :file:`vacumm-scripts/tutorials` subdirectory in the
      installed package directory (see :meth:`get_lib_dir`).
    - Else in the :file:`scripts/tutorials`
      subdirectory of the main distribution tree
      (see :meth:`get_dist_dir`).

    .. warning:: It raises an :class:`VACUMMError` error if not found
        and if ``raiseerr`` is True.
    """
    return get_scripts_dir('tutorials', raiseerr=raiseerr)


def get_home_conf_dir():
    """Get the directory which contains the user config dir

    It defaults to :file:`~/.config` and may be set with the :envvar:`XDG_CONFIG_HOME`
    environmement variable.
    """
    base_dir = os.environ.get('XDG_CONFIG_HOME')
    if not base_dir:
        base_dir = os.path.expanduser(os.path.join('~', '.config'))
    return base_dir

def get_user_conf_dir():
    """Get directory that contains user configuration files

    This directory contains main and secondary configuration files.
    These files are usually referenced in the general configuration file.
    On linux, it should be :file:`~/.config/vacumm`.
    It may be set with the :envvar:`VACUMM_CONFIG_DIR`.
    """
    conf_dir = os.environ.get('VACUMM_CONFIG_DIR')
    if conf_dir is None:
        conf_dir = os.path.join(get_home_conf_dir(), 'vacumm')
    return conf_dir

def get_user_conf_file(fname=None):
    """Get main user configuration file

    :Params:

        - **fname**, optional: Relative or absolute file name.
          It may be set with the :envvar:`VACUMM_CONFIG_FILE`.
          The config directory is prepended if the file name
          is relative.

    Shortcut for::

        os.path.join(get_user_conf_dir(), fname)

    :See also: :func:`get_user_conf_dir`, :func:`edit_user_conf_file`
    """
    if fname is None: fname = os.environ.get('VACUMM_CONFIG_FILE')
    if fname is None: fname = 'vacumm.cfg'
    if os.path.isabs(fname): return fname
    return os.path.join(get_user_conf_dir(), fname)

def get_com_conf_dir():
    """Get directory of common alternate configuration files
    These files are not included in the distribution.
    They have a similar role to user config files but they
    are accessible to all users.
    This directory can be at two different places, depending on if the library is
    an installed version or developers version.
    - If installed : in :file:`vacumm-config` subdirectory in the
      installed package directory (see :meth:`get_lib_dir`).
    - Else in the :file:`config` subdirectory of the main distribution tree
      (see :meth:`get_dist_dir`).
    """
    # Installed librairy
    lib_dir = get_lib_dir()
    conf_dir = os.path.join(lib_dir, 'vacumm-config')
    if os.path.exists(conf_dir):
        return conf_dir

    # Distributed library (dev)
    dist_dir = get_dist_dir()
    if dist_dir is not None:
        return os.path.join(dist_dir, 'config')

    # Fall back to the use config dir
    return get_user_conf_dir()

def get_com_conf_file():
    """Get common user configuration file
    Shortcut for::
        os.path.join(get_com_conf_dir(), 'vacumm.cfg')
    :See also: :func:`get_com_conf_dir`, :func:`get_user_conf_file`
    """
    return os.path.join(get_user_conf_dir(), 'vacumm.cfg')


def get_dir_dict():
    """Get the following directory names as a dictionary

    - ``lib_dir`` (see :func:`get_lib_dir`)
    - ``data_dir`` (see :func:`get_data_dir`)
    - ``dist_dir`` (see :func:`get_dist_dir`)
    - ``tut_dir`` (see :func:`get_tut_dir`)
    - ``scripts_dir`` (see :func:`get_scripts_dir`)
    - ``user_conf_dir`` (see :func:`get_user_conf_dir`)
    - ``com_conf_dir`` (see :func:`get_com_conf_dir`)
    """

    return dict(data_dir=get_data_dir(), lib_dir=get_lib_dir(),
        tut_dir=get_tut_dir(), dist_dir=get_dist_dir(),
        com_conf_dir=get_com_conf_dir(),
        user_conf_dir=get_user_conf_dir(), scripts_dir=get_scripts_dir())

def get_dl_dir(suggest=None, filedesc=None, quiet=False):
    """Get the download directory"""
    # Try what is suggested
    if suggest is not None:
        if not os.access(suggest, os.W_OK|os.R_OK|os.X_OK):
            try:
                os.makedirs(suggest)
                return suggest
            except:
                pass
        else:
            return suggest

    # Already chosen?
    dl_dir = vacumm.VACUMM_CFG['vacumm'].get('dl_dir')
    if dl_dir and os.path.exists(dl_dir): return dl_dir

    # Download
    user_data_dir = os.path.join(get_user_conf_dir(), 'data')
    if quiet: # Into config dir
        dl_dir = user_data_dir
    else: # Ask

        # Possible directories
        home_data_dir = os.path.expanduser(os.path.join('~', 'data'))
        dl_dirs = [user_data_dir, home_data_dir]
        vacumm_data_dir = get_data_dir()
        if os.access(vacumm_data_dir, os.W_OK|os.R_OK|os.X_OK):
            dl_dirs.append(vacumm_data_dir)
        nc = max([len(p) for p in dl_dirs])

        # Message
        if filedesc is None: filedesc = 'missing data file'
        msg = 'Please specifiy or choose a directory where to download %s:\n'%filedesc
        msg += '- specify it -> enter a valid path\n'
        for i, p in enumerate(dl_dirs):
            default = '*' if i==0 else ''
            msg += '- %s -> %i%s\n'%(p.ljust(nc), i+1, default)
        msg += "? "

        # Ask choice
        choice = None
        while True: #in range(1, len(dl_dirs)+1):
            choice = input(msg).strip()
            if len(choice)==0:
                choice = '1'
                print('1')
            if len(choice)==1 and choice.isdigit(): # choose within the list of dirs
                choice = int(choice)
                if choice<1 or choice>nc:
                    print("Wrong choice: %i"%choice)
                    continue
                break
            elif choice: # specifiy dir
                if not os.path.isdir(choice): # does not exist
                    try:
                        os.makedirs(choice)
                    except:
                        print("Can't create directory: "+choice)
                        continue
                if not os.access(choice, os.W_OK|os.R_OK|os.X_OK): # no access
                    print("No write access to directory: "+choice)
                    continue
                break
        if isinstance(choice, int):
            dl_dir = dl_dirs[choice-1]
        else:
            dl_dir = choice

    # Remember it
    vacumm.VACUMM_CFG['vacumm']['dl_dir'] = dl_dir
    rem = input('Remember this choice? [Y/n] ').strip()
    if not rem.lower().startswith('n'):
        save_config_value(vacumm.VACUMM_CFG['vacumm'], 'dl_dir', dl_dir)

    return dl_dir


#: Config specifications file
VACUMM_CFGSPECS_FILE = os.path.join(os.path.dirname(__file__),  'vacumm.ini')

#: Config specifications
VACUMM_CFGSPECS = ConfigObj(VACUMM_CFGSPECS_FILE, interpolation=False, list_values=False)


class VCValidator(Validator):

    def __init__(self, *args, **kwargs):

        self._vc_checkmode = kwargs.pop('checkmode', 'warn')
        self._vc_warnmsg = kwargs.pop('warnmsg',
                'Bad configuration value: {value}. '
                'Error: "{error}". '
                'It must conform the follwing specs: {check}. '
                'Switching to default value: {default}.')
        Validator.__init__(self, *args, **kwargs)

    def check(self, check, value, missing=False, warnmsg=None):
        """
        Unlike the original validator, it replaces the value by its default
        when an error occurs.
        """
        fun_name, fun_args, fun_kwargs, default = self._parse_with_caching(check)

        if missing:
            if default is None:
                # no information needed here - to be handled by caller
                raise VdtMissingValue()
            value = self._handle_none(default)

        if value is None:
            return None

        if self._vc_checkmode=='normal':
            return self._check_value(value, fun_name, fun_args, fun_kwargs)
        try:
            return self._check_value(value, fun_name, fun_args, fun_kwargs)
        except Exception as  e:
            if default is None:
                raise VdtMissingValue()
            error = e.message
            if warnmsg is None:
                warnmsg = self._vc_warnmsg
            warn(warnmsg.format(**locals()))
            return self._check_value(default, fun_name, fun_args, fun_kwargs)

#: VCValidator default instance
VACUMM_VDT = VCValidator(checkmode='warn',
    warnmsg='Bad configuration value: {value}. '
        'Error message: "{error}". '
        'It must conform the follwing specs: {check}. '
        'Switching to default value: {default}. '
        'Please refer to the documentation to fix your configuration: '
        'http://wwww.ifremer.fr/vacumm/user.install.config.html')

#: VCValidator special instance to check individual values
VACUMM_VDT_LIVE = VCValidator(checkmode='warn',
    warnmsg='Bad value: {value}. '
        'Error message: "{error}". '
        'It must conform the follwing specs: {check}. '
)

def load_cfg(cfgfile=None, merge=True, live=False, validate=True):
    """Load the configuration

    :Params:

        - **cfgfile**, optional: A cfg file. Defaults to the result of
          :func:`get_user_conf_file`. Set it to False to not load
          the default user config file (:func:`get_user_conf_file`).
        - **merge**, optional: Merge it with the currently loaded
          configuration :data:`VACUMM_CFG`.
        - **live**, optional: Load but does not store it into store in
            :data:`VACUMM_CFG`. So it must used from the return argument.
        - **validate**, optional: Validate the configuration using
          :class:`VCValidator`. It checks that options conform specifications
          and fallbacks to default value when an error is detected or when
          the option is missing.

    :Return: A :class:`configobj.ConfigObj` instance.
    """

    # Load
    if cfgfile=='user' or cfgfile is None:
        cfgfile = get_user_conf_file()
    elif cfgfile=='com':
        cfgfile = get_com_conf_file()
    elif cfgfile=='default':
        cfgfile = False
    if cfgfile is None:
        cfgfile = get_user_conf_file()
    if (not isinstance(cfgfile, ConfigObj)  or not isinstance(cfgfile, str) or
            not os.path.exists(cfgfile)):
        cfgfile = ""
#        warn('Invalid cfgfile passed to load_cfg. Skipping.')
    cfg = ConfigObj(cfgfile, configspec=VACUMM_CFGSPECS, interpolation='template')

    # Default section
    cfg.setdefault('DEFAULT', get_dir_dict())

    # Validate
    if validate:
        cfg.validate(VACUMM_VDT)

    # Merge with currently loaded config?
    if not hasattr(vacumm, 'CFG') or vacumm.CFG is None: # nothing to merge with
        merge = False
    if merge:
        oldcfg = vacumm.CFG
        if live:
            oldcfg.copy()
        oldcfg.merge(cfg)
        cfg = oldcfg

    # Store it?
    if not live:
        vacumm.VACUMM_CFG = vacumm.CFG = cfg

    return cfg

def _get_parent_sections_(sec):
    secs = []
    while sec.parent is not sec:
        secs.append(sec.name)
    return secs[::-1]

def _get_sections_(sec):
    """Get current section and parent names

    :Params:

        - **sec**: A :class:`~configobj.Section` instance or a single section name
          or a list of section names.

    :Return: ``sec, sections``

        - ``sec``: Current :class:`~configobj.Section`.
        - ``sections``: Parent section names.
    """
    if isinstance(sec,  str):
        sec = [sec]
    if isinstance(sec, list):
        sections = sec
        sec = vacumm.VACUMM_CFG
        for secname in sections:
            sec = sec[secname]
    else:
        sections = _get_parent_sections_(sec)
    return sec, sections

def get_config_value(sec, option):
    """Get a config value

    :Params:

        - **sec**: A :class:`~configobj.Section` instance or a single section name
          or a list of section names.
        - **option**: Option name.

    :Return: ``value``
    """
    sec, sections = _get_sections_(sec)
    return sec[option]

def set_config_value(sec, option, value):
    """Set a config value

    .. note:: This function set the value of un optin for current session usage
        only. If you want this option value to be persistent, please use
        function :func:`save_config_value`.

    :Params:

        - **sec**: A :class:`~configobj.Section` instance or a single section name
          or a list of section names.
        - **option**: Option name.
        - **value**: Option value.

    :Return: ``value``
    """
    sec, sections = _get_sections_(sec)
    sec[option] = value
    return value


def save_config_value(sec, option, value):
    """Save a single config option to the user config file

    It load the user config file without validation, change the option,
    and write the config back to the user file.
    """
    # Save the live value
    sec[option] = value

    # Load and update the user config
    cfg = load_cfg(live=True, validate=False)
    sec = cfg
    for secname in _get_parent_sections_()[:-1]:
        sec = sec[secname] = {}
    sec[option] = value

    # Save it
    cfg.write(get_user_conf_file())


def get_cfg_checked(sec, option, value=None):
    """Check a value or get its default value"""
    # Current section and parent names
    sec, sections = _get_sections_(sec)

    # Default value
    if value is None:
        return sec[option]

    # Find the check specs and check!
    sec = VACUMM_CFGSPECS
    for secname in sections:
        sec = sec[secname]
    check = sec[option]
    return VACUMM_VDT_LIVE.check(value, check)

def get_cfg_path(sec, option, expand=True):
    """Get and interpret a config value as a path with expansions"""
    if not option in sec:
        return
    path = sec[option]
    if expand:
        for op in (os.path.expanduser, os.path.expandvars):
            path = op(path)
    return path

def _print_header_(text, nc):
    print('#'*nc)
    print('##  {:nc-4}'.format(text))
    print('#'*nc)

def print_config(section=None, system=True, direc=True, config=True, packages=True,
        extended=False, user=True, headers=True):
    """Print current configuration

    :Params:

        - **system**: print the base of the config.
        - **direc**: print key directories of vacumm.
        - **config**: print the vacumm config.
        - **packages**: print the version of some important packages.
        - **extended**: force printing everything.

    This is equivalent to (if ``extended=False``)::

        >>> print get_config(section, user, asstring=True)

    """
    nc = 70
    if extended: 
        config = packages = True
    if user:
        if headers: _print_header_('System', nc)
        print('Python executable: '+sys.executable)
        print('Python version: '+' '.join(sys.version.splitlines()))
        print('Platform: '+sys.platform)
    if direc:
        if headers: _print_header_('VACUMM directories', nc)
        for dd in get_dir_dict().items():
            print('%s: %s'%dd)
    if config:
        if headers: _print_header_('VACUMM configuration', nc)
        cfg = vacumm.VACUMM_CFG
        if section and section in vacumm.VACUMM_CFG:
            cfg = vacumm.VACUMM_CFG[section]
        print(cfg)
    if packages:
        if headers: _print_header_('Versions', nc)
        from .__init__ import __version__
        print('VACUMM: '+__version__)
        for mname in 'numpy', 'matplotlib', ('mpl_toolkits.basemap', 'basemap'), 'scipy':
            try:
                if isinstance(mname, tuple):
                    m = getattr(__import__(mname[0]), mname[1])
                else:
                    m = __import__(mname)
                v = m.__version__
            except ImportError as e:
                v = 'ERROR: '+e.message
            if isinstance(mname, tuple): mname = mname[0]
            print('%s: %s'%(mname.split('.')[-1], v))


def get_default_config():
    """
    Shortcut to::

        load_cfg(cfgfile=False, merge=False, live=False)

    """
    return load_cfg(cfgfile=get_com_conf_file(), merge=False, live=True)

def write_default_config(cfgfile):
    """Write the defaults config (see :func:`get_default_config`) to a file"""
    get_default_config().write(cfgfile)


def check_data_file(section,
        quiet=None, suffix=None, avail=False, check_license=True):
    """Check the existence of a file whose path is stored in a configuration file

    The path is given by ``section["file"]``.
    If the file is not found, it may download it using an url
    given by ``section["url"]``

    :Tasks:

        #. If ``suffix`` is a list, call itself in a loop with ``suffix``
           set to each element of the list.
        #. Get the ``path``, ``url`` and ``license`` of the data file.
        #. If ``avail is True``, only check the existence and return the status.
        #. If the file is not found (and ``quiet is False``), ask the user
           for it. If the path specified is empty, simply go further.
           If is not empty, update the configuration and return
           it, else raise an error.
        #. If a license info is present, tell the user that he must
           have the authorization to download it. If the user does not agree,
           raise an error.
        #. Ask the user where to download the file using :func:`get_dl_dir`.
        #. Download the file using :func:`download_file`.
        #. Update the configuration if needed.
        #. Return the path.

    :Params:

        - **section**: A section object with "file", "url" and "license" options.
        - **quiet**, optional: Don't ask question about where to download the
          file, don't display any info, don't ask for authorization.
        - **suffix**, optional: A suffix or a list of suffixes to append to
          the path and url before using it.
        - **check_license**, optional: If a license info is found,
          show ot and ask the user if he has autorization for downloading it.
        - **avail**, optional: Check availability only (see below).

    :Return:

        - A single path or a list of paths.
        - ``None`` or a list of them if nothing found.
        - If ``avail=True``:

            - ``0``: Data file not on disk and not downloadable
              (however, user can specify the path manually).
            - ``1``: File is not on disk, but downloadable.
            - ``2``: File is on disk.

    :Examples:

        >>> check_data_file(
            VACUMM_CFG['vacumm.bathy.shorelines']['shapefiles']['histolitt'],
            suffix=['.shp', '.dbf'])
    """
    if quiet is None:
        quiet = not sys.stdin.isatty()

    # Loop on suffixes
    if isinstance(suffix, (list, tuple)):
        paths = []
        for i, suf in enumerate(suffix):
            paths.append(check_data_file(section, quiet=quiet, suffix=suf,
                check_license=i==0))
        return paths
    elif not isinstance(suffix, six.string_types):
        suffix = ''

    # Get local and remote file names
    path = section['file']
    url = section['url']
    license = section['license']
    path += suffix
    if url is not None: url += suffix

    # Only check availability
    if os.path.exists(path):
        if avail: return 2
        return path
    if avail:
       return 0 if url is None else 1

    # Ask for this path
    if not quiet:
        print("Can't find data file: "+path)
        guessed = input("If you know what is this file and where it is on your system,\n"
            "please enter its name here, or leave it empty: \n").strip()
        if guessed:
            if not os.path.exists(guessed):
                print('File not found')
            else:
                save_config_value(section, 'file', path[:len(path)-len(suffix)])
                return guessed

    # Check url
    if not url:
        if quiet: return
        raise VACUMMError('Data file not found and not url provided for downloading it: '+path)

    # License for downloading
    if license and check_license:
        nc = len(license)+4
        lic = '#'*nc+'\n'
        lic += '# %s #\n'%license
        lic += '#'*nc
        print("VACUMM is about to download this data file: %s\n"%url + \
            "We suppose you have requested the authorization and are not responsible for your choice.\n" + \
            "If you're not sure you are allowed to do it, please abort.\n" + \
            "For more information about this data file and associated distribution license, check this:\n"+lic)
        while True:
            try:
                c = input("Would you like to download it? [y/N]\n")
            except:
                c = 'n'
            if not c: c='n'
            if c.startswith('n'):
                raise VACUMMError("Download interrupted -> can't access to data")
            if c.startswith('y'): break

    # Download directory
    dl_dir = path and os.path.dirname(path)
    if not path or not os.access(dl_dir, os.W_OK|os.R_OK|os.X_OK):
        dl_dir = get_dl_dir(quiet=quiet, suggest=dl_dir)

    # Download
    basename = url.split('/')[-1]
    dl_path = os.path.join(dl_dir, basename)
    download_file(url, dl_path, quiet=quiet)

    # Fix configuration
    if path!=dl_path:
        save_config_value(section, 'file', dl_path[:len(path)-len(suffix)])
    return dl_path



# UTILITIES

class URLLibProgressBar:
    """Progress bar compatible with :func:`urllib.urlretrieve`

    It looks like::

        [=======>        22%                  ]

    :Example:

        >>> import urllib
        >>> mybar = URLProgress(width=50)
        >>> urllib.urlretrieve("http://web.net/file", "myfile", reporthook=mybar)

    :Inspired from: http://code.activestate.com/recipes/168639-progress-bar-class/
    """

    def __init__(self, width=80):
        self.prog = "[]"
        self.width = max(width, 6)
        self.amount = -1
        self.np = -1
        self(0, 1, 1)  # Init

    def __call__(self, count, blocksize, totalsize):
        amount = int(count*blocksize*100./totalsize)
        if amount==self.amount: return
        end =  amount>=100
        amount = min(amount, 100)
        pwidth = self.width - 2
        np = max(int(round((amount / 100.0) * pwidth)), 1)
        if self.np==np: return
        self.np = np
        if np == pwidth:
            prog = "[%s]" % ('='*pwidth)
        else:
            prog = "[%s>%s]" % ('='*(np-1), ' '*(pwidth-np))
        percentPlace = (len(prog) / 2) - len(str(amount))
        percentString = str(amount) + "%"
        prog = ''.join([prog[0:percentPlace], percentString,
            prog[percentPlace+len(percentString):]])
        sys.stdout.write("\r%s"%prog)
        if end: sys.stdout.write("\n")



def download_file(url, locfile, quiet=False):
    """Download a file with a progress bar using :func:`urllib.urlretrieve`"""
    locfile = os.path.realpath(locfile)
    if not quiet:
        print('Downloading %s to %s'%(url, locfile))
    dl_dir = os.path.dirname(locfile)
    if not os.path.exists(dl_dir): os.makedirs(dl_dir)
    six.moves.urllib.request.urlretrieve(url, locfile, reporthook=URLLibProgressBar() if not quiet else None)
    return locfile


def edit_file(fname, editor=None, verbose=True):
    """Edit a file

    :Params:

        - **fname**: File name.
        - **editor**, optional: Editor to use (default to :envvar:`VISUAL` or
          :envvar:`EDITOR`.

    :Source: Mercurial code
    """
    # File name
    if not os.path.exists(fname):
        raise IOError('File not found: '+fname)

    # Guess editor
    if editor is None:
        editor = os.environ.get("VISUAL") or os.environ.get("EDITOR", "vi")

    # Run editor
    if verbose: print('Editing user configuration file: '+fname)
    if subprocess.call(shlex.split("%s \"%s\"" % (editor, fname))):
        raise VACUMMError("Error editing file: "+fname)
    if verbose: print('End of file edition')

def edit_user_conf_file(fname='vacumm.cfg', editor=None, verbose=True):
    """Edit the user configuration file

    :Params:

        - **fname**, optional: Relative file name, which defaults to :file:`vacumm.cfg`.
        - **editor**, optional: Editor to use (default to :envvar:`VISUAL` or
          :envvar:`EDITOR`.

    Examples:

        >>> edit_user_conf_file()
        >>> edit_user_conf_file('my_other_config.cfg', editor='gedit')

    :See also: :func:`edit_file` and :func:`get_user_conf_file`
    """
    edit_file(get_user_conf_file(fname), editor=editor, verbose=verbose)


def handle_print_args(args):
    """Handle argparse args from ``vacumm_config.py print ...``"""

    # Options
    if not (args.system or args.direc or args.config or args.packages):
        args.system = args.direc = args.config = args.packages = True

    # Print
    print_config(system=args.system, direc=args.direc, config=args.config,
        packages=args.packages, headers=not args.no_header)

def handle_edit_args(args):
    """Handle argparse args from ``vacumm_config.py edit ...``"""

    edit_user_conf_file(editor=args.editor)





