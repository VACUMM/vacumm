# -*- coding: utf8 -*-
"""Utilities to manage VACUMM basic configuration of modules"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
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


import os, sys
import re
import glob
import subprocess,  shlex
import numpy.distutils
from ConfigParser import SafeConfigParser, ConfigParser
import StringIO
from vacumm import VACUMMError

def get_lib_dir():
    """Directory of the :mod:`vacumm` package"""
    return os.path.abspath(os.path.dirname(__file__))


def get_dist_dir():
    """Upper directory of VACUMM distribution tree or ``None`` for an installed package"""
    lib_dir = get_lib_dir()
    dist_dir = os.path.abspath(os.path.join(lib_dir, '../../..'))
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
    if subdir and isinstance(subdir, basestring):
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

#def get_conf_dir():
#    """Get directory of secondary configuration files
#
#    This directory contains secondary configuration files.
#    These files are usually referenced in the general configuration file.
#    It can be at two different places, depending on if the library is
#    an installed version or developers version.
#
#    - If installed : in :file:`vacumm-config` subdirectory in the
#      installed package directory (see :meth:`get_lib_dir`).
#    - Else in the :file:`config` subdirectory of the main distribution tree
#      (see :meth:`get_dist_dir`).
#    """
#    # Installed librairy
#    lib_dir = get_lib_dir()
#    conf_dir = os.path.join(lib_dir, 'vacumm-config')
#    if os.path.exists(conf_dir):
#        return conf_dir
#
#    # Distributed library (dev)
#    dist_dir = get_dist_dir()
#    if dist_dir is not None:
#        return os.path.join(dist_dir, 'config')
#
#    # Fall back to the use config dir
#    return get_user_conf_dir()

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

def get_com_conf_file():
    """Get common user configuration file


    Shortcut for::

        os.path.join(get_com_conf_dir(), 'vacumm.cfg')

    :See also: :func:`get_com_conf_dir`, :func:`get_user_conf_file`
    """
    return os.path.join(get_user_conf_dir(), 'vacumm.cfg')

def get_dir_dict(modname=None):
    """Get the following directory names as a dictionary

    - ``lib_dir`` (see :func:`get_lib_dir`)
    - ``data_dir`` (see :func:`get_data_dir`)
    - ``dist_dir`` (see :func:`get_dist_dir`)
    - ``tut_dir`` (see :func:`get_tut_dir`)
    - ``com_conf_dir`` (see :func:`get_com_conf_dir`)
    - ``user_conf_dir`` (see :func:`get_user_conf_dir`)
    - ``scripts_dir`` (see :func:`get_scripts_dir`)
    - ``mod_dir`` (directory of ``modname``)

    :Params:

        - **modname**, optional: Module name.
    """

    dd = dict(data_dir=get_data_dir(), lib_dir=get_lib_dir(),
        tut_dir=get_tut_dir(), dist_dir=get_dist_dir(),
        com_conf_dir=get_com_conf_dir(), conf_dir=get_user_conf_dir(),
        user_conf_dir=get_user_conf_dir(), scripts_dir=get_scripts_dir())
    if modname is not None:
        if isinstance(modname, basestring):
            if os.path.isdir(modname):
                mod_dir = modname
                modname = None
            else:
                try:
                    exec 'import %s'+modname
                    modname = eval(modname)
                except:
                    mod_dir = os.path.join(get_lib_dir(), *modname.split('.')[1:-1])
                    modname = None
        if modname is not None:
            mod_dir = os.path.dirname(modname.__file__)
        dd.update(mod_dir=mod_dir)
    return dd

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
    dl_dir = get_config_value('vacumm', 'dl_dir')
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
            choice = raw_input(msg).strip()
            if len(choice)==0:
                choice = '1'
                print '1'
            if len(choice)==1 and choice.isdigit(): # choose within the list of dirs
                choice = int(choice)
                if choice<1 or choice>nc:
                    print "Wrong choice: %i"%choice
                    continue
                break
            elif choice: # specifiy dir
                if not os.path.isdir(choice): # does not exist
                    try:
                        os.makedirs(choice)
                    except:
                        print "Can't create directory: "+choice
                        continue
                if not os.access(choice, os.W_OK|os.R_OK|os.X_OK): # no access
                    print "No write access to directory: "+choice
                    continue
                break
        if isinstance(choice, int):
            dl_dir = dl_dirs[choice-1]
        else:
            dl_dir = choice

    # Remember it
    rem = raw_input('Remember this choice? [Y/n] ').strip()
    if not rem.lower().startswith('n'):
        set_config_value('vacumm', 'dl_dir', dl_dir)

    return dl_dir

def get_config_files(section=None, check=False, user=1):
    """Get the list of configuration files that can be loaded

    It lists all configuration files that should be loaded
    to search for a configuration value. It first list
    default configuration files, then add user specific files.


    :Params:

        - **section**, optional: Target section of the configuration files.
          It is generally a module name.
          In this case, user and default configuration files
          located in the where the module is are inserted
          into the list.
          It can also be a list of module name.
          Eventually, if the special value ``"__all__"`` is passed, it searches
          for all configuration files named :file:`vacumm.cfg`
          in the whole library tree.
        - **user**, optional: If ``False`` or ``0``, load only default config.
          If ``1`` or ``True``, append user configuration files
          to the list. If ``2``, list only user config files.


    :Example:

        Suppose you are in directory :file:`/home/me/myproject`,
        :func:`get_config_files` is called from module
        :mod:`vacumm.bathy.shorelines`, and your use VACUMM
        as a developper, i.e. the library is not installed but used
        directly in the installation tree whose path is
        :file:`/home/me/python/vacumm`:

        >>> get_config_files('vacumm.bathy.shorelines')

        The following list is returned :

            #. /home/me/python/vacumm/lib/python/vacumm/bathy/shorelines/vacumm.cfg (defaults)
            #. /home/me/python/vacumm/lib/python/vacumm/vacumm-config/vacumm.cfg (defaults)
            #. /home/me/.config/vacumm/vacumm.cfg

        **Content of last files overides content of first files.**

        Other usages:

        >>> get_config_files('__all__', user=False)
        >>> get_config_files(['toto.tutu', 'tutu.tata'])
    """
    # Some directories
    lib_dir = get_lib_dir()
    com_conf_dir = get_com_conf_dir()
    user_conf_dir = get_user_conf_dir()
    mod_dirs = []
    if section=='__all__':
        for root, dirs, files in os.walk(lib_dir):
            if '.svn' in dirs:
                dirs.remove('.svn')
            for d in dirs:
                if d.startswith('vacumm-') or \
                        not os.path.exists(os.path.join(root, d, '__init__.py')):
                    dirs.remove(d)
            if 'vacumm.cfg' in files:
                mod_dirs.append(root)
    elif section:
        if hasattr(section, '__name__'):
            section = section.__name__
        if not isinstance(section, list):
            section = [section]
        for sec in section:
            if sec.startswith('vacumm.') or sec=='vacumm':
                subnames = sec.split('.')[1:] # skip "vacumm"
                if len(subnames)>1:
                    subnames = subnames[:-1] # skip file name
                mod_dirs.append(os.path.join(lib_dir, *subnames))
    dist_dir = get_dist_dir()

    # Configuration files to scan (no check)
    cfgfiles = []
    if user!=2:
        if mod_dirs:
            cfgfiles.extend([os.path.join(mod_dir, 'vacumm.cfg') for mod_dir in mod_dirs])  # module level / default
    if user:
        cfgfiles.append(os.path.join(com_conf_dir, 'vacumm.cfg')) # global tree / user
        cfgfiles.append(os.path.join(user_conf_dir, 'vacumm.cfg')) # global tree / user

    # check existence
    if check:
        cfgfiles = [f for f in cfgfiles if os.path.exists(f)]
    return cfgfiles

def get_config_value(section, option, regexp=False, ispath=False, user=True,
    expand=True, vars=None, raw=False, umode='auto', cfg=None, getcfg=False,
    parent_section=None, parent_option=None):
    """Get the a config value in available VACUMM configuration files

    The list of possible configuration files is provided by :func:`get_config_files`.

    Configuration values can use expanded variables in two ways:

        1. Internal option name are automatically expanded to their value
           if found in the form ``%(option_name)s``. In addition to internal
           options, you can provide your own (var_name,value) pairs by
           passing a dictionary to the ``vars`` keyword.
           Some external options pointing to directory names are provided
           automatically (but can be overidden):

            - ``lib_dir`` (see :func:`get_lib_dir`)
            - ``data_dir`` (see :func:`get_data_dir`)
            - ``dist_dir`` (see :func:`get_dist_dir`)
            - ``tut_dir`` (see :func:`get_tut_dir`)
            - ``user_conf_dir`` or ``conf_dir`` (see :func:`get_user_conf_dir`)
            - ``scripts_dir`` (see :func:`get_scripts_dir`)
            - ``mod_dir`` (directory of the module when ``section`` is a module name)

        2. User home ``~`` and environement variable (in the form ``${HOME}``
           are expanded if keyword ``expand`` is set to ``True``.


    :Params:

        - **section**: A section in the config, typically a full module name
          (see :func:`get_config_files`)
        - **option**: Option to read in the section.
        - **regexp**, optional: ``option`` is interpreted a regular expression.
          If ``True``, it returns a dictionary of (option,value) items found.
        - **ispath**, optional: if ``True``, value is converted using
          :func:`numpy.distutils.misc_util.allpath`.
        - **vars**, optional: Dictionary used to expand variables in the form
          ``%(var_name)s``.
        - **user**, optional: Also scan user files. If ``raw`` is True, value
          with ``conf_dir`` expanded using ``user_conf_dir``, and if differents values
          are found, they are treated depending on ``mode``.
        - **umode**, optional: If ``user`` is True, and two possible values are found
          (one with ``conf_dir``, one with ``user_conf_dir``), they are treated
          depending on ``umode``:

            - ``"merge"``: Value found using ``user_conf_dir`` has priority
              over value found using ``conf_dir``.
            - ``"all"``: A two-element list is returned if elements are different.
            - Else, if value seems to be a config file (ending with ".ini"
              or ".cfg"), ``mode="all"`` else ``mode="merge"``.
        - **cfg**, optional: External :class:`ConfigParser.ConfigParser` instance
          or an external config file. If not provided, configuration is given
          by :func:`get_config` (internal configuration).

    :Returns: If ``regexp`` is `` False``, return the value or None if not found,
        else return a dictionary.

    :Example:

        >>> path = get_config_value('vacumm.bathy.shorelines', 'shapefile_histolitt')
        >>> paths = get_config_value('vacumm.bathy.shorelines', 'shapefile_.*', regexp=True)
        >>> paths = get_config_value('my.module', vars=dict(mydir='/home/me/dir'))
    """
    ispath=False # buggy

    # Load configuration(s)
    if isinstance(cfg, (basestring, list)):
        cfgfile = cfg
        cfg = SafeConfigParser()
        cfg.read(cfgfile)
    elif not isinstance(cfg, ConfigParser):
        cfg = get_config(section=section, user=user, parent_section=parent_section,
            parent_option=parent_option)

    # Always have a tuple of configurations
    cfgs = (cfg, ) if not isinstance(cfg, tuple) else cfg

    # Directory names (for expansion)
    kwget = get_dir_dict(section)
    kwget['conf_dir'] = '(%s|%s)'%(kwget['com_conf_dir'], kwget['user_conf_dir'])
    if vars:
        kwget.update(vars)

    # Get value(s)
    if not regexp:

        # Collect
        values = [(cfg.has_option(section, option) and
                cfg.get(section, option, vars=kwget, raw=raw) or None) for cfg in cfgs]

        # Detect tuple of values
        values = _opt2tuple_(values)

        # Filter
        if user:
            value = _filter_uvalues_(values, umode)
        else:
            value = values[0]

    else: # Only options that matches regular expression

        # Collect
        value = {}
        for cfg in cfgs:
            if cfg.has_section(section):
                re_match = re.compile(option, re.I).match
                for opt in cfg.options(section):
                    m = re_match(opt)
                    if m:
                        value.setdefault(opt, [])
                        value[opt].append(cfg.get(section, opt, vars=kwget, raw=raw))

        # Filters
        for opt, values in value.items():

            # Detect tuple of values
            values = _opt2tuple_(values)

            if user:
                value[opt] = _filter_uvalues_(values, umode)
            else:
                value[opt] = values[0]

    # Operator on values
    if value:

        # Operators
        ops = []
        # - fix path
        if ispath:
            ops.append(numpy.distutils.misc_util.allpath)
        # - expand user name and environment variables
        if expand:
            ops.extend([os.path.expanduser, os.path.expandvars])

        # Apply
        if ops:
            value = _ops_on_values_(value, *ops)

    return (value, cfgs) if getcfg else value

def _ops_on_values_(values, *ops):
    """Apply operators on tuple, list, dict or scalar values"""

    # Dict
    if isinstance(values, dict):
        for key, val in values.items():
            values[key] = _ops_on_values_(values[key], *ops)
        return values

    # Scalar
    if not isinstance(values, (tuple, list)):
        for op in ops:
            values = op(values)
        return values

    # Tuple, list
    return type(values)(_ops_on_values_(v, *ops) for v in values)

def _filter_uvalues_(values, umode):
    """Filter values depending on umode

    umode is guess if "auto"

    values is a N-element list:

        - first: default value
        - others: user values
    """
    # Auto mode = 'all' if config files as values
    if umode!='all' and umode!='merge':
        umode = 'merge'
        for value in values:
            if value is not None and _re_cfgfile(value):
                umode = "all"
                break

    # Remove Nones
    values = filter(None, values)

    # Merge values (take last)
    if umode=="merge":
        return None if not values else values[-1]

    # Unique
    vv = []
    for v in values:
        if v not in vv:
            vv.append(v)

    return tuple(vv)

def get_config_sections(cfg=None, parent_section=None, parent_option=None, getcfg=False):
    """Get sections names

    :Params:

        - **parent_section/option**, optional: Section (module name) and option of
          the main vacumm configuration that give the path to the target configuration
          file. If not provided, it uses `cfg` the main vacumm configuration.
        - **getcfg**, optional: Also return the configurations object.
    """
    cfgs = get_config(cfg=cfg, parent_section=parent_section, parent_option=parent_option)
    if not isinstance(cfgs, tuple): cfgs = cfgs,
    sections = []
    for cfg in cfgs:
        sections.extend(cfg.sections())
    if getcfg: return sections, cfg
    return sections


def load_configs(cfgfiles):
    """Read configuration files in merged and separate configuations

    :Params:

        - **cfgfiles**: Configuration files.

            - If a file, simply return the associated :class:`~ConfigParser.SafeConfigParser`
              instance.
            - If a list of files, does the same merging all the files.
            - If a tuple, create a configuration for each of its elements.

    :Example:

        >>> cfg = load_configs('myfile.cfg')
        >>> cfg = load_configs(['myfile1.cfg', 'myfile2.cfg'])
        >>> cfg1, cfg2 = load_configs(('vacumm.cfg', ['myfile1.cfg', 'myfile2.cfg']))

    """
    single = not isinstance(cfgfiles, tuple)
    if single:
        cfgfiles = (cfgfiles, )
    cfgs = ()
    for cfgfile in cfgfiles:
        if isinstance(cfgfile, basestring):
            cfgfile = [cfgfile]
        cfg = SafeConfigParser()
        cfg.read(cfgfile)
        cfgs += cfg,
    if single: return cfgs[0]
    return cfgs




def get_config(section=None, user=True, asstring=False, cfg=None,
    parent_section=None, parent_option=None, split=False):
    """Get a :class:`ConfigParser.SafeConfigParser` instance of current configuration or as a string

    :Params:

        - **section**, optional: A module or its name, or ``"__all__"``
          to get a full config (config files in module drectories are also scanned).
        - **user**, optional: Also scan available user configuration files.
        - **asstring**, optional: Get it as string.
        - **parent_section/option**, optional: Section (module name) and option of
          the main vacumm configuration that give the path to the target configuration
          file. If not provided, it uses `cfg` the main vacumm configuration.

    :See also: :func:`get_config_files`
    """
    # Already there?
    if not isinstance(cfg, ConfigParser): # No

        # List of files to scan
        if isinstance(cfg, basestring): # Single file name
            cfgfiles = [cfg]

        elif isinstance(cfg, (list, tuple)): # A list of files
            cfgfiles = cfg

        elif parent_section and parent_option: # File name from parent config
            cfgfiles = get_config_value(parent_section, parent_option, user=user, umode='all')
            if isinstance(cfgfiles, list):
                cfgfiles = tuple(cfgfiles)

        else: # General config files
            cfgfiles = get_config_files(section=section, user=False)
            if user:
                cfgfiles = (cfgfiles, get_com_conf_file(), get_user_conf_file())

        # Load configurations
        cfg = load_configs(cfgfiles)


    # Return config
    if not asstring:
        return cfg

    # As string
    if isinstance(cfg, tuple):
        cfg = merge_configs(cfg)
    S = StringIO.StringIO()
    cfg.write(S)
    s = S.getvalue()
    S.close()
    return s

_re_opt2tuple_match = re.compile(r"^(?:\(([^|)]+)\|([^|)]+)(?:\|([^|)]+))?(?:\|([^|)]+))?\))([^)]*)$").match
def _opt2tuple_(values):
    """Split "(a|b)c" into ("ab","bc")"""
    if values is None: return
    if not isinstance(values, (list, tuple)):
        values = [values]
    out = []
    for value in values:
        if value is None: continue
        m = _re_opt2tuple_match(value)
        if m is not None:
            g = list(m.groups())
            vv = filter(None, g[:-1])
            if g[-1] is not None:
                for i, v in enumerate(vv):
                    vv[i] += g[-1]
            out.extend(vv)
        else:
            out.append(value)
    return tuple(out)


_re_cfgfile = re.compile(r'\.(cfg|ini)$').search
def merge_configs(cfgs, forprint=True, raw=True):
    """Merge configurations

    Options of the last config overrides options of the other ones,
    except for config file names (values then end with .ini or .cfg)
    that are store as a tuple.

    .. warning::

        This makes output config incompatible with
        :meth:`~ConfigParser.ConfigParser.getint`, etc.
    """
    # Init
    if not isinstance(cfgs, (list, tuple)):
        cfgs = [cfgs]
    cfg = SafeConfigParser()

    # Loop on configs
    kwget = get_dir_dict()
    kwget['mod_dir'] = '%%(mod_dir)s'
    pending = {}
    for cfgi in cfgs:
#        print '-'*50
#        S = StringIO.StringIO()
#        cfgi.write(S)
#        s = S.getvalue()
#        S.close()
#        print s

        if not isinstance(cfgi, SafeConfigParser):
            cfgfiles = cfgi
            cfgi = SafeConfigParser()
            cfgi.read(cfgfiles)

        cfg._defaults.update(cfgi.defaults())

        for section in cfgi.sections():

            if section not in cfg.sections():
                cfg.add_section(section)
            pending.setdefault(section, {})

            for option, rvali in cfgi.items(section, raw=raw):

                if option in pending[section]: # already identified as an option of config files
                    if rvali not in pending[section][option]:
                        pending[section][option] = tuple(pending[section][option]) + (rvali,)
                elif _re_cfgfile(cfgi.get(section, option, vars=kwget)):
                    pending[section][option] = rvali,
                else:
                    cfg.set(section, option, rvali)

    # Convert pending options
    for section, options in pending.items():
        for option, value in options.items():
            if forprint:
                if len(value)==1:
                    value = value[0]
                else:
                    value = '(%s)'%', '.join(value)
            else:
                value = '(%s)'%'|'.join(value)
#            print section, option, value
            cfg.set(section, option, value)

    return cfg

def _print_header_(text, nc):
    print '#'*nc
    print '##  %s'%text.ljust(nc-4)
    print '#'*nc

def print_config(section='__all__', system=True, direc=True, config=True, packages=True,
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
    if extended: base = config = packages = True
    if user:
        if headers: _print_header_('System', nc)
        print 'Python executable: '+sys.executable
        print 'Python version: '+' '.join(sys.version.splitlines())
        print 'Platform: '+sys.platform
    if direc:
        if headers: _print_header_('VACUMM directories', nc)
        for dd in get_dir_dict().items():
            print '%s: %s'%dd
    if config:
        if headers: _print_header_('VACUMM configuration', nc)
        print get_config(section=section, user=user, asstring=True).strip()
    if packages:
        if headers: _print_header_('Versions', nc)
        from __init__ import __version__
        print 'VACUMM: '+__version__
        for mname in 'numpy', 'matplotlib', ('mpl_toolkits.basemap', 'basemap'), 'scipy':
            try:
                if isinstance(mname, tuple):
                    m = getattr(__import__(mname[0]), mname[1])
                else:
                    m = __import__(mname)
                v = m.__version__
            except ImportError, e:
                v = 'ERROR: '+e.message
            if isinstance(mname, tuple): mname = mname[0]
            print '%s: %s'%(mname.split('.')[-1], v)


def get_default_config(section='__all__'):
    """
    Shortcut to::

        get_config(section='__all__', user=False)

    """
    return get_config(section=section, user=False)

def write_default_config(cfgfile):
    """Write the defaults config (see :func:`get_default_config`) to a file"""
    f = open(cfgfile, 'w')
    get_default_config().write(f)
    f.close()

def set_config_value(section, option, value=None, quiet=False, cfgfile=None):
    """Create or update user configuration file

    The directory of this file is given by :func:`get_user_conf_dir`
    and its name is :file:`vacumm.cfg`. Therefore it should be:
    :file:`~/.config/vacumm/vacumm.cfg`.

    :Params:

        - **section**: A module or its name.
        - **option**: Option name.
        - **value**: Value of the option. If ``None``, option is removed.
        - **cfgfile**, optional: Configuration file. If not provided,
          internal user configuration file is used
          (:file:`vacumm.cfg` in user config directory -- given by :func:`get_user_conf_dir`).

    :Example:

        >>> set_config_value('vacumm.bathy.bathy', 'cfgfile_gridded',
            '%(mod_dir)s/bathy.gridded.cfg') # set
        >>> import vacumm.bathy.bathy as B
        >>> set_config_value(B, 'cfgfile_gridded') # remove

    :Output: Path of the configuration file

    """
    if hasattr(section, '__name__'):
        section = section.__name__

    # Load or check
    if cfgfile is None:
        cfgfile = os.path.join(get_user_conf_dir(), 'vacumm.cfg')
    if not os.access(cfgfile, os.W_OK):
        if os.path.exists(cfgfile):
            raise VACUMMError("Can't write to config file: "+cfgfile)
        else:
            udir = os.path.dirname(cfgfile)
            try:
                os.makedirs(udir)
            except:
                raise VACUMMError("Can't write to config file: "+cfgfile)
    cfg = SafeConfigParser()
    cfg.read(cfgfile)

    # Update
    if value is not None:
        value = str(value)
        if not cfg.has_section(section):
            cfg.add_section(section)
        cfg.set(section, option, value)
        if not quiet:
            print 'Updated user configuration file (%s) with:'%cfgfile
            print ' [%(section)s]\n  %(option)s=%(value)s'%locals()
    elif cfg.has_option(section, option):
        cfg.remove_option(section, option)
        if not quiet:
            print 'Removed the following option from user configuration file (%s):'%cfgfile
            print ' [%(section)s]\n  %(option)s'%locals()

    # Save
    f = open(cfgfile, 'w')
    cfg.write(f)
    f.close()
    return cfgfile

def check_data_file(section, option, parent_section=None, parent_option=None,
    quiet=None, suffix=None, avail=False, check_license=True):
    """Check the existence of a file whose path is stored in a configuration file

    Two cases are treated:

        1. Path value is accessible from the vacumm configuration (:func:`get_config`),
           using ``section`` (module name) and ``option``.
        2. Path value is accessible from a secondary configuration file using
           ``parent_section`` (module name) and ``parent_option``, and the
           path value of this config file is accessible from the vacumm configuration
           using ``section`` and ``name``.

    If the file is not found, it may download it using an url whose value is
    accessible at the same place as the path and with the same option
    name + '_url'.

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

        - **section**: Section where to find the path to data.
        - **option**: Option to read the path.
        - **parent_section**, optional: Section of the vacumm configuration
          where to find the path to the secondary configuration file that has
          the path stored within.
        - **parent_option**, optional: Option of the vacumm configuration
          to read the path to the secondary configuration.
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

        >>> check_data_file('vacumm.bathy.shorelines', 'shapefile_histolitt', suffix=['.shp', '.dbf'])
        >>> check_data_file('etopo2', 'file', parent_section='vacumm.bathy.bathy',
        ... parent_option='cfgfile_gridded')
    """
    if quiet is None:
        quiet = not sys.stdin.isatty()
    # Loop on suffixes
    if isinstance(suffix, (list, tuple)):
        paths = []
        for i, suf in enumerate(suffix):
            paths.append(check_data_file(section, option, parent_section=parent_section,
                parent_option=parent_option, quiet=quiet, suffix=suf, check_license=i==0))
        return paths
    elif not isinstance(suffix, basestring):
        suffix = ''

    # Get local and remote file names
    path, cfg = get_config_value(section, option, umode='merge', getcfg=True,
        parent_section=parent_section, parent_option=parent_option)
    url = get_config_value(section, option+'_url', umode='merge', cfg=cfg)
    license = get_config_value(section, option+'_license', umode='merge', cfg=cfg)
    if path is None:
        raise VACUMMError("Can't determine path of data file. Here are the config specs:\n"
            "section='%(section)s',  option='%(option)s', parent_section='%(parent_section)s', parent_option='%(parent_option)s'"%locals())
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
        print "Can't find data file: "+path
        guessed = raw_input("If you know what is this file and where it is on your system,\n"
            "please enter its name here, or leave it empty: \n").strip()
        if guessed:
            if not os.path.exists(guessed):
                print 'File not found'
            else:
                _set_config_path_(guessed, section, option, parent_section,
                    parent_option, parent and cfgpath, suffix)
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
        print "VACUMM is about to download this data file: %s\n"%url + \
            "We suppose you have requested the authorization and are not responsible for your choice.\n" + \
            "If you're not sure you are allowed to do it, please abort.\n" + \
            "For more information about this data file and associated distribution license, check this:\n"+lic
        while True:
            try:
                c = raw_input("Would you like to download it? [y/N]\n")
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
        _set_config_path_(dl_path, section, option, parent_section,
            parent_option, parent_section and cfgpath, suffix)
    return dl_path

def _set_config_path_(path, section, option, parent_section, parent_option, cfgpath, suffix):
    """Save path into config"""
    path = path[:len(path)-len(suffix)]
    if parent_option and parent_section.startswith('vacumm.'):
        if isinstance(cfgpath, list): cfgpath = cfgpath[-1]
        if True: #not os.access(cfgpath,  os.W_OK|os.R_OK): # new config file
            cfgpath = os.path.join(get_user_conf_dir(),
                os.path.basename(cfgpath))
        set_config_value(section, option, path, cfgfile=cfgpath) # save data path
        set_config_value(parent_section, parent_option, cfgpath) # save config path
    elif section.startswith('vacumm.'):
        set_config_value(section, option, path)
    return path


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

import urllib, sys

def download_file(url, locfile, quiet=False):
    """Download a file with a progress bar using :func:`urllib.urlretrieve`"""
    locfile = os.path.realpath(locfile)
    if not quiet:
        print 'Downloading %s to %s'%(url, locfile)
    dl_dir = os.path.dirname(locfile)
    if not os.path.exists(dl_dir): os.makedirs(dl_dir)
    urllib.urlretrieve(url, locfile, reporthook=URLLibProgressBar() if not quiet else None)
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
    if verbose: print 'Editing user configuration file: '+fname
    if subprocess.call(shlex.split("%s \"%s\"" % (editor, fname))):
        raise VACUMMError("Error editing file: "+fname)
    if verbose: print 'End of file edition'

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


if __name__=='__main__':
#    print_config()


    import vacumm.config2;reload( vacumm.config2); print vacumm.config2.get_config_value('vacumm.bathy.bathy','cfgfile_gridded')

