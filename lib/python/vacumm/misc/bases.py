#!/usr/bin/env python
# -*- coding: utf8 -*-
#
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

__author__ = 'Jonathan Wilkins'
__email__ = 'wilkins@actimar.fr'
__doc__ = '''

This module is intented to provide features like logging, configuration and
debugging into a base class.

The main tool of this module is the `class:`Object` class which
may be used in high level components of your applications by inheriting from this
base class.


.. warning::

    You must be aware that all `class:`Object` instances hold
    a logger and an configuration manager.
    Therefore you should not inherit from `class:`Object` in all of
    your specialized classes especially when your application instanciate and holds
    a lot of instances of this class.
    The `class:`Object` class should only be used in main classes
    of your application.



'''
 
import inspect, os, pdb, pprint, sys, types

import configobj

from vacumm import vacumm_warning
from vacumm.misc.config import ConfigManager
from vacumm.misc.exception import getDetailedExceptionInfo
from vacumm.misc.misc import kwfilter, dict_merge
from vacumm.misc.log import Logger, logger, get_str_levels
from argparse import ArgumentParser
from optparse import OptionParser

# Test application origin: a python script or an executable (py2exe/cx_freeze)
isfreezed = os.path.basename(sys.executable) == os.path.basename(sys.argv[0])

# Program identity variables
prog = os.path.realpath(isfreezed and sys.executable or sys.argv[0])
if prog in ['', '-', '-c']:
    progname = prog = 'python'
    progdir = os.getcwd()
else:
    progname = os.path.basename(prog)
    progdir = os.path.dirname(prog)


def get_prog_name(noext=False):
    '''
    Get the currently running program name

    :Params:
        - **noext**: If True, remove extension from program name

    :Return:
        - str: The program name, or "python" in interactive interpreter mode.
    '''
    return noext and os.path.splitext(progname)[0] or progname


def get_prog_dir():
    '''Get the currently running program's directory.'''
    return progdir


def func_name(iframe=0):
    '''
    Get the name of the calling function

    Parameters:
        - **iframe**: int: Get the iframe'th caller function name.
    '''
    #return inspect.currentframe().f_back.f_code.co_name
    return sys._getframe(1+iframe).f_code.co_name


def _set_ext_(fname, ext):
    """Change a file extension

    Parameters:
        - **fname**: File path.
        - **ext**: Remove extension if False, replace it if string,
          or leave it.
    """
    if ext is False or isinstance(ext, basestring):
        fname, suf = os.path.splitext(fname)
        if ext is not False:
            if isinstance(ext, basestring):
                suf = ext
            if not '.' in suf:
                suf = '.'+suf
            fname = fname + suf
    return fname

def code_file_name(iframe=0,  ext=True):
    '''
    Get the name of the file that hosts the code where it is called

    Parameters:
        - **iframe**: int: Get the iframe'th caller function name.
        - **ext**: Remove extension if False, replace it if string,
          or leave it.
    '''
    fname = sys._getframe(1+iframe).f_code.co_filename
    fname = _set_ext_(fname, ext)
    return fname

def code_base_name(iframe=0, ext=True):
    '''
    Get the basename of the file that hosts the code where it is called

    Parameters:
        - **iframe**: int: Get the iframe'th caller function name.
        - **ext**: Remove extension if False, replace it if string,
          or leave it.
    '''
    return os.path.basename(code_file_name(iframe+1, ext))

def code_dir_name(iframe=0):
    '''
    Get the dirname of the file that hosts the code where it is called

    Parameters:
        - **iframe**: int: Get the iframe'th caller function name.
        - **ext**: Remove extension if False, replace it if string,
          or leave it.
    '''
    return os.path.dirname(code_file_name(iframe+1))

def stack_trace(iframe=0):
    def etrace(frame, pad=60):
        o = frame.f_locals.get('self', '') # not the perfect way...
        if o: o = '%s(%s).'%(o.__class__.__name__, id(o))
        else:
            o = frame.f_locals.get('cls', '') # not the perfect way...
            if o: o = '%s.'%(o.__name__)
        o += frame.f_code.co_name
        return ('%%-%ds'%(pad)+'  (%s:%s)')%(o, frame.f_code.co_filename, frame.f_lineno)
    iframe = 1 + iframe
    stack = inspect.stack()
    frame = stack[iframe][0]
    traces = []
    last = ' '+etrace(frame)
    while True:
        try:
            iframe += 1
            frame = stack[iframe][0]
            traces.append(' %s'%(etrace(frame)))
        except: break
    traces.reverse()
    traces.append(last)
    return 'Stack trace (innermost last):\n'+'\n'.join(traces)


try:
    import psutil
    _psutil_process = psutil.Process(os.getpid())
    def psinfo():
        '''Get cpu and memory usage string (require **psutil** module)'''
        if hasattr(_psutil_process, 'get_cpu_percent'):
            cpu = _psutil_process.get_cpu_percent()
            mem = _psutil_process.get_memory_percent()
            meminfo  = _psutil_process.get_memory_info()
        else:
            cpu = _psutil_process.cpu_percent()
            mem = _psutil_process.memory_percent()
            meminfo = _psutil_process.memory_info()
        rss, vsz = meminfo[:2]
        ps = '[CPU: %d%%  MEM: %d%%  RSS: %dMo  VSZ: %dMo'%(cpu, mem, rss / 2**20, vsz / 2**20)
        if len(meminfo) >= 6:
            ps += '  SHR: %dMo  DAT: %dMo'%(meminfo[2] / 2**20, meminfo[5] / 2**20)
        ps += ']'
        return ps
except Exception, e:
    Logger.default.verbose('psinfo disabled: %s', e)
    psinfo = lambda:'Ressources informations not available (no module psutil)'


def describe(obj, stats=None, format=pprint.pformat):
    '''Return the object decription depending on its type.

    Usefull with numpy and cdms variables/axes

    :Params:
        - **obj**: The object to describe.
        - **stats**: If True, include numerical information like min, max, mean, count, ...

    :Return:
        - The object's summary string
    '''
    try:
        import cdms2, MV2, numpy
        from cdms2.avariable import AbstractVariable
        from cdms2.axis import AbstractAxis
        if isinstance(obj, configobj.ConfigObj):
            obj = obj.dict()
        if not isinstance(obj, (cdms2.avariable.AbstractVariable, cdms2.axis.AbstractAxis, numpy.ndarray)):
            return format(obj)
        otype = obj.__class__.__name__ # type(obj)
        sh, sz = MV2.shape(obj), MV2.size(obj)
        mi = ma = av = co = None
        if stats and sz:
            try:
                if hasattr(obj, 'typecode') and obj.typecode() not in ('c',):
                    mi, ma, av, co = MV2.min(obj), MV2.max(obj), MV2.average(obj), MV2.count(obj)
            except: logger.exception('Error getting statistics of object %s', type(obj))
        if isinstance(obj, AbstractVariable):
            return '%s: %s, shape: (%s), order: %s%s'%(
                otype, obj.id,
                ','.join('%s=%s'%(a.id, a.shape[0]) for a in obj.getAxisList()),
                obj.getOrder(),
                stats and ', min: %s, max: %s, avg: %s, count: %s'%(mi, ma, av, co) or '')
        elif isinstance(obj, AbstractAxis):
            return '%s: %s, size: %s%s'%(
                otype, obj.id, sz,
                stats and ', min: %s, max: %s, avg: %s, count: %s'%(mi, ma, av, co) or '')
        else:#if isinstance(obj, numpy.ndarray):
            return '%s: shape: %s%s'%(
                otype, sh,
                stats and ', min: %s, max: %s, avg: %s, count: %s'%(mi, ma, av, co) or '')
    except Exception, e:
        logger.exception('Error getting description of object %s', type(obj))
        return '%s (error getting description)'%(type(obj))


################################################################################
# http://stackoverflow.com/questions/5892619/static-and-instance-methods-in-python
################################################################################

class _classinstancemethod_wrapper(object):
    def __init__(self, func, obj, type):
        self.func = func
        self.obj = obj
        self.type = type
    def __call__(self, *args, **kw):
        assert not kw.has_key('self') and not kw.has_key('cls'), (
            "You cannot use 'self' or 'cls' arguments to a "
            "classinstancemethod")
        return self.func(*((self.obj, self.type) + args), **kw)
    def __repr__(self):
        if self.obj is None:
            return ('<bound class method %s.%s>'
                    % (self.type.__name__, self.func.func_name))
        else:
            return ('<bound method %s.%s of %r>'
                    % (self.type.__name__, self.func.func_name, self.obj))

class classinstancemethod(object):
    """
    Decorator which acts like a class method when called from a class, like an
    instance method when called by an instance.  The method should
    take two arguments, 'self' and 'cls'; one of these will be None
    depending on how the method was called.
    """
    def __init__(self, func):
        self.func = func
    def __get__(self, obj, type=None):
        return _classinstancemethod_wrapper(self.func, obj=obj, type=type)

################################################################################


class Class(type):
    '''
    Vacumm's base metaclass used to initialize implementing classes.

    The implementing classes (in fact via the :class:`Object` class below) can have methods
    called at the class definition (at present, the :meth:`Object._init_class()` classmethod).

    Usefull to initialize class attributes such like class configuration, logging, debugging...
    '''
    # More details in Object class comments below

    # <fix>
    # A fix for (at least) doc generation because the inspect module throws
    # an AttributeError in a getattr call to get this attribute...
    # This attribute seems to be used by the abc module (abstract base classes)
    # I do not met some side effects using this fix yet.
    __abstractmethods__ = []
    # </fix>

    # Not yet needed
    #def __new__(mcs, name, bases, dct):
    #    return type.__new__(mcs, name, bases, dct)

    def __init__(cls, name, bases, dct):
        type.__init__(cls, name, bases, dct)
        cls._init_class(name, bases, dct)


_logging_proxies = list(map(lambda f: (f.lower(),f.lower()), get_str_levels()+['exception']))
_logging_proxies += [
    ('get_loglevel', 'get_level_name'),
    ('set_loglevel', 'set_level'),
    ('is_verbose', 'is_verbose'),
    ('is_debug', 'is_debug'),
]
def add_logging_proxies(cls):
    '''
    Register some Logger shortcuts for a given class (cls.<method> => cls.get_logger().<method>):

        - %s

    These shortcuts will be bound to the class and its instances using its get_logger method
    which must return a :class:`Logger` instance, so this method must be callable from the class or its
    instances (you may use the :class:`classinstancemethod` decorator as of :meth:`Object.get_logger`)
    '''
    def wrap_logging_function(cls , cfunc, lfunc=None):#, addextra=False):
        if not lfunc: lfunc = cfunc
        @classinstancemethod
        def wrapper(obj, cls, *a, **k):
#            if addextra:
#                k = k.copy()
#                if lfunc != 'exception:'
#                    k.setdefault('extra', {})
#                    # We cannot overwrite funcName so add a new keyword
#                    k['extra'].setdefault('funcname', func_name(1))
            return getattr(obj.get_logger() if obj is not None else cls.get_logger(), lfunc)(*a,**k)
            #return getattr(obj._logger if obj is not None else cls._class_logger, lfunc)(*a,**k)
        wrapper.__name__ = cfunc
        wrapper.__doc__ = '''Wrapper to the %(lfunc)r method of this object logger (see :meth:`Logger.%(lfunc)s` )'''%vars()
        setattr(cls, cfunc, wrapper)
        
        # Register code to be skipped when logging
        cls.get_logger().skipCaller(wrapper.func)
        # XXX because _classinstancemethod_wrapper is created on classinstancemethod.__get__(),
        # we directly use _classinstancemethod_wrapper.__call__
        cls.get_logger().skipCaller(_classinstancemethod_wrapper.__call__.__func__.__code__)
        
    for args in _logging_proxies:
        wrap_logging_function(cls, *args)
add_logging_proxies.__doc__ %= ('\n        - '.join(map(lambda f: '%s => %s'%f, _logging_proxies)),)

class Object(object):
    '''
    Vacumm's base class proving common usefull features:

        - configuration management
            - keeping trace of parent class configs
            - allowing per-instance config usage (default to class config)
        - logging
        - debugging
            - tracking (pdb)
            - exceptions details
            - process memory usage

    Most features work at class level to reduce variables in instances.

    This class is designed to be implemented by the main working classes,
    but there is no real restriction about that.

    The default configuration is automatically loaded when importing the class module, this
    behavior is also applied for subclasses in other modules.

    Instance initial configuration (:func:`get_config`) is a copy of the class
    configuration (:func:`get_default_config`) made at instance creation time.

    '''
    # All subclasses will use this metaclass
    __metaclass__ = Class
    # For example, the following classes definition:
    #   In a vacumm.a.py module:
    #     from vacumm import Object
    #     class A(Object): pass
    #   In b vacumm.b.py module:
    #     import a
    #     class B(a.A): pass
    #
    # will produce the following methods calls at class definition (import):
    #
    #   Class.__new__ (in Class): cls: <class 'vacumm.Class'>
    #   Object.__init__ (in Class): self: <class 'vacumm.Object'>
    #   Class.__new__ (in Class): cls: <class 'vacumm.Class'>
    #   A.__init__ (in Class): self: <class 'a.A'>
    #   Class.__new__ (in Class): cls: <class 'vacumm.Class'>
    #   B.__init__ (in Class): self: <class 'b.B'>
    #
    # and then at instanciation (b=B()):
    #
    #   B.__new__ (in Object): cls: <class 'b.B'>
    #   B.__init__ (in Object): cls: <b.B object at 0x2ad3b50>
    #   A.__init__
    #   B.__init__
    #

    _log_level = 'info'
    _cfg_debug = False
    _log_obj_stats = False

    # Not yet needed
    #def __new__(cls, *args, **kwargs):
    #    #print '%s.%s (in Object): cls: %s, args: %s, kwargs: %s'%(cls.__name__, func_name(), cls, args, kwargs)
    #    #print '%s.%s (in Object): cls: %s'%(cls.__name__, func_name(), cls)
    #    #
    #    return object.__new__(cls, *args, **kwargs)

    # ==========================================================================
    # Class initializer (not instance initializer !!!)
    # ==========================================================================
    @classmethod
    def _init_class(cls, name, bases, dct):
        '''
        Class initialization method, called when the class is defined.
        '''
        # Setup class logging
        cls._class_logger = Logger(name=name, name_filters=[name])
        add_logging_proxies(cls)
        # Configuration management
        if not hasattr(cls, '__config_managers'):
            cls.__config_managers = {}
        if not hasattr(cls, '__config_defaults'):
            cls.__config_defaults = {}
        cls.load_default_config(nested=True)
        cls.init_class(name, bases, dct)

    @classmethod
    def init_class(cls, name, bases, dct):
        '''Redefine this method if you need class initialization'''
        pass

    # ==========================================================================
    # Logging features
    # ==========================================================================

    @classmethod
    def get_class_logger(cls):
        'Return the :class:`Logger` instance bound to this class.'
        return cls._class_logger

    @classmethod
    def set_class_logger(cls, logger):
        'Set the :class:`Logger` instance bound to this class.'
        cls._class_logger = logger

    @classinstancemethod
    def get_logger(obj, cls):
        '''
        Return the :class:`Logger` instance bound to this class or instance,
        according to the way this method is called (from class or instance).
        '''
        return cls._class_logger if obj is None else obj._logger

    @classinstancemethod
    def set_logger(obj, cls, logger):
        '''
        Set the :class:`Logger` instance bound to this class or instance,
        according to the way this method is called (from class or instance).
        '''
        if obj is None: cls._class_logger = logger
        else: obj._logger = logger

    logger = property(lambda o,*a,**k:o.get_logger(*a,**k), lambda o,*a,**k:o.set_logger(*a,**k), None,
        'A :class:`Logger` instance. You can use this object for all logging '
        'operations related to this class')

    # ==========================================================================
    # Configuration features
    # ==========================================================================

    @classmethod
    def get_config_spec_file(cls):
        '''Return (and define) the class specification file path'''
        try:
            cfgfile = '%s.ini'%(os.path.splitext(os.path.realpath(inspect.getfile(cls)))[0])
        except TypeError: # occure if creating a subclass in interactive python shell
            cfgfile = None
        return cfgfile

    @classmethod
    def get_parent_config_spec(cls):
        '''Get the merged config specifications of all parents'''
        cfg = None
        for c in cls.__bases__:
            if not hasattr(c, 'get_config_spec'): continue
            cs = c.get_config_spec()
            if cfg is None:
                cfg = cs
            else:
                cfg = dict_merge(cfg, cs)
        return cfg


    @classmethod
    def get_config_spec(cls):
        '''Load the config specs as ConfigObj object

        It merges the specs of the current class and those of parents classes
        '''
        cfg = None
        spec = cls.get_config_spec_file()

        # If specification (and so defaults) file defined and exists
        if spec and os.path.isfile(spec):

            # Load (temporary) the file
            cfg = configobj.ConfigObj(spec, list_values=False, interpolation=False)

            # NOTE: list_values=False, interpolation=False are set because list values
            #       are handled by the manager (parse error otherwise)

            # If a config section lookup is defined and present, load the section
            sec = cls.get_config_section_name()
            if sec and sec in cfg:
                #cfg = configobj.ConfigObj(cfgspec[sec], list_values=False, interpolation=False)
                cfg = configobj.ConfigObj(cfg[sec])

        # Merge with parents
        pcfg = cls.get_parent_config_spec()
        if cfg is None:
            cfg = pcfg
        elif pcfg is not None:
            cfg = dict_merge(cfg, pcfg, mergesubdicts=False)
        return cfg


    @classmethod
    def get_config_section_name(cls):
        '''Return (and define) the class specification section name'''
        return cls.__name__

    @classmethod
    def get_config_manager(cls, reload=False, encoding=None):
        '''
        Get the configuration manager for this class (cls).
        This manager and its underlying configuration specification
        must not be dynamically changed as it is fixed at design time.

        .. note:: this method is also the config manager lazy loader

        '''
        # Populate this class config manager if not yet done or reload request
        if not cls in cls.__config_managers or reload:

            # Get the ConfigObj object of specifications
            cfg = cls.get_config_spec()

            # NOTE: If no spec / no class section, class use empty spec
            cfgmgr = ConfigManager(cfg, encoding=encoding)
            cls.__config_managers[cls] = cfgmgr
            if cls._cfg_debug:
                cls.debug('Loaded config manager of class %s with spec:\n  %s', cls.__name__, '\n  '.join(cfgmgr._configspec.write()))

        return cls.__config_managers[cls]

    @classmethod
    def load_default_config(cls, config=None, nested=None, apply=True, encoding=None):
        '''Load / update the class (unique) default configuration'''
        cfgmgr = cls.get_config_manager(encoding=encoding)
        cfgdef = cfgmgr.defaults()
        cfgsec = cls.get_config_section_name()
        if isinstance(nested, basestring):
            cfgsec = nested
        # Change the class default config if required
        if config is not None:
            cfg = configobj.ConfigObj(config, encoding=encoding)
            # If a config section lookup is required
            # Otherwise, the whole passed config will be taken
            if nested:
                if cfgsec in cfg:
                    cfg = configobj.ConfigObj(cfg[cfgsec], encoding=encoding)
                # Section not found, use empty config
                else:
                    cfg = configobj.ConfigObj(encoding=encoding)
            cfgdef = cfgmgr.load(cfg)
        cls.__config_defaults[cls] = cfgdef
        if apply:
            cls.apply_default_config()
        if cls._cfg_debug:
            cls.debug(
                'Loaded %s default configuration:'
                '\n  section:  %s'
                '\n  nested:   %s'
                '\n  from:     %s'
                '\n  loaded:   '
                '\n  %s',
                cls, cfgsec, nested, config, '\n  '.join(cls.get_default_config().write()))

    @classmethod
    def get_default_config(cls, encoding=None):
        '''Get the default configuration (copy)'''
        return configobj.ConfigObj(cls.__config_defaults[cls], encoding=encoding)

    @classmethod
    def get_default_config_str(cls, encoding=None):
        '''Get the default configuration as a string'''
        return '\n  '.join(configobj.ConfigObj(cls.__config_defaults[cls], encoding=encoding).write())

    @classmethod
    def apply_default_config(cls, config=None, encoding=None):
        '''
        This will turn on debug mode of various features (config, objects stats)
        if the loglevel of nested Logger configuration is debug.

        Set the default log level for the newly created objects based on nested
        Logger configuration.

        Subclasses may override this to apply/update according to the new config

        :Params:
            - **config**: The new config loaded by :meth:`load_default_config()`.

        .. note::
            - overriding this method will obviously shunt its default beahvior, you'll then have to call original method if needed

        '''
        if config is None: config = cls.get_default_config(encoding=encoding)
        #cls.get_logger().load_config(config, nested=True)
        loglvl = logger.get_level_name() #loglvl = cls.get_logger().get_level_name().lower()
        isdbg = logger.is_debug() #isdbg = loglvl == 'debug'
        cls._log_level = loglvl
        cls._cfg_debug = config.get('cfg_debug', isdbg or cls._cfg_debug)
        cls._log_obj_stats = config.get('log_obj_stats', isdbg or cls._log_obj_stats)

    @classmethod
    def from_config(cls, config, *args, **kwargs):
        '''Create a cls instance using args and kwargs and load config.

        The **nested** named argument (in kwargs) is extracted before
        creating the instance and then passed to load_config.

        :Params:
            - **config**: A configuration file (str) or object (ConfigObj).
            - **args** and **kwargs**: Passed to the object constructor, without parmeters described above.

        :Return:
            - The created object of class cls

        '''
        loadkw = dict(((a,kwargs.pop(a)) for a in ('nested', 'apply') if a in kwargs))
        obj = cls(*args, **kwargs)
        obj.load_config(config, **loadkw)
        return obj

    def load_config(self, config=None, nested=None, apply=True, cfgpatch=None, encoding=None, **kwargs):
        '''Load / update the instance configuration

        :Params:
            - **config**: A configuration file (str) or object (ConfigObj) or
              None to load defaults, or an :class:`~argparse.ArgumentParser` object.
            - **nested**: Load from a nested config section instead of the whole config.
                          If True, use the section name returned by :meth:`get_config_section_name()`
                          Else if a string, use the section name defined by the **nested** string
            - **cfgpatch**: A manual patch to apply to the config once loaded.
            - Other options are passed to
              :meth:`~vacumm.misc.config.ConfigManager.arg_parse`.

        :Return: A :class:`ConfigObj` object or :class:`ConfigObj`,options tuple if
            an :class:`~argparse.ArgumentParser` object has been passed
        '''
        mgr = self.get_config_manager(encoding=encoding)
        sec = self.get_config_section_name()
        if not hasattr(self, '_config'):
            self._config = self.get_default_config(encoding=encoding)
        self._options = None
        if config is not None:
            if isinstance(nested, basestring):
                    sec = nested
            if isinstance(config, ArgumentParser):
                self._config, self._options = mgr.arg_parse(config,
                    nested=nested and sec or nested, getargs=True, **kwargs)
            else:
                cfg = configobj.ConfigObj(config, interpolation=False, encoding=encoding)
                # If a nested section lookup is required
                # Otherwise, the whole passed config will be taken
                if nested and sec and sec in cfg: # If not found, self._config remain unchanged
                    cfg = configobj.ConfigObj(cfg[sec], interpolation=False, encoding=encoding)
                self._config = mgr.load(cfg)
        if cfgpatch is not None:
            if not isinstance(cfgpatch, list):
                cfgpatch = [cfgpatch]
            for patch in cfgpatch:
                mgr.cfg_patch(self._config, patch)
        if apply:
            self.apply_config(self._config)
        if self._cfg_debug:
            self.debug(
                'Loaded %s configuration:'
                '\n  section:  %s'
                '\n  nested:   %s'
                '\n  from:     %s'
                '\n  loaded:   '
                '\n  %s',
                self.__class__, sec, nested, config, '\n  '.join(self._config.write()))
        return self._config

    def save_config(self, outfile=None, nested=None):
        if isinstance(outfile, basestring):
            outfile, close = file(outfile, 'w'), True
        else:
            close = False
        config = configobj.ConfigObj(self._config)
        if nested:
            if isinstance(nested, basestring):
                sec = nested 
            else:
                sec = self.get_config_section_name()
            # XXX config.dict() is required, otherwise section will not be correctly written ([] missing)
            config = configobj.ConfigObj({sec:config.dict()})
        r = config.write(outfile)
        if close:
            outfile.close()
        return r

    def get_options(self):
        """Get :attr:`options`"""
        return getattr(self, '_options', None)
    options = property(fget=get_options, doc='Options loaded from the commandline parser or None')

    def get_config(self, copy=True):
        '''Get the instance's config'''
        if copy:
            return configobj.ConfigObj(self._config)
        return self._config
    config = property(fget=get_config, doc='Current configuration')

    def get_config_str(self):
        '''Get the instance's config as a string'''
        return '\n  '.join(configobj.ConfigObj(self._config).write())

    def apply_config(self, config):
        '''Subclasses may override this to apply/update according to the new config

        :Params:
            - **config**: The new config loaded by :meth:`load_config()`.

        .. note::
            - overriding this method will obviously shunt its default beahvior,
              you'll then have to call original method if needed

        '''
        self.logger.load_config(config, nested=True)
        # ...


    # ==========================================================================
    # Debugging/Introspection features
    # ==========================================================================

    @classmethod
    def func_name(cls, iframe=0):
        __doc__ = func_name.__doc__
        return func_name(iframe+1)

    @staticmethod
    def stack_trace(iframe=0):
        __doc__ = stack_trace.__doc__
        return stack_trace(iframe+1)

    @classmethod
    def exception_trace(cls):
        '''Return a huge detailed exception traceback'''
        return getDetailedExceptionInfo()

    @classmethod
    def trace(cls, iframe=0, iftty=True):
        '''
        Start pdb debugger

        :Params:
            - **iframe**: frame index entry point of the debugger, relative to the caller
            - **iftty**: if True, disable this call in a non interactive execution

        .. note::
            - For debugging purpose only: do not let trace calls in a production
              environment, even if an interactive test is done !
        '''
        if sys.stdin.isatty():# and sys.stdout.isatty():
            pdb.Pdb().set_trace(sys._getframe(1+iframe))

    @classmethod
    def describe(cls, obj, **kwargs):
        kw = dict(stats=cls._log_obj_stats, format=cls.pformat)
        kw.update(kwargs)
        return describe(obj, **kw)

    @classmethod
    def pformat(cls, obj, indent=2, width=80, depth=None):
        '''Pretty print an object'''
        # default is: indent=1, width=80, depth=None
        return pprint.pformat(obj, indent, width, depth)

    # ==========================================================================
    # Other stuff
    # ==========================================================================

    @classmethod
    def kwfilter(cls, kwargs, filters=None, *args, **kwa):
        '''
        Shortcut to :func:`vacumm.misc.misc.kwfilter` with a filters argument which
        defaults to this class lowercase name.
        '''
        if filters is None: filters = cls.__name__.lower()
        return kwfilter(kwargs, filters, *args, **kwa)

    # ==========================================================================
    # Initializer
    # ==========================================================================
    def __init__(self, **kwargs):
        '''

        :Keyword arguments:
            - **config*: load a configuration (:meth:`load_config()`)
            - **logger**: A :class:`Logger` instance.
            - **logger_<param>**: ``<param>`` is passed to the
              :class:`Logger` constructor. Not used if **logger** is passed.

        '''
        # Setup logging
        lkw = kwfilter(kwargs, 'logger', {'name':self.__class__.__name__})
        if isinstance(lkw.get('config', None), Object):
            lkw['config'] = lkw['config'].get_logger()
        lkw['name_filters'] = list(lkw.get('name_filters', [])) + [self.__class__.__name__]
        if 'logger' in lkw:
            self._logger = lkw['logger']
            if not isinstance(self._logger, Logger):
                vacumm_warning(self.__class__.__name__.split('.')[-1]+
                    '{} is initialised with an invalid logger type')
        else:
            self._logger = Logger(**lkw)
        self._logger.skipCaller(self.get_class_logger().skipCaller())
        # Load passed or default configuration
        self.load_config(kwargs.get('config', None),  cfgpatch=kwargs.get('cfgpatch', None))


################################################################################
# http://code.activestate.com/recipes/204197-solving-the-metaclass-conflict/
################################################################################

import inspect, types, __builtin__

###################### preliminary: two utility functions ######################

def skip_redundant(iterable, skipset=None):
    'Redundant items are repeated items or items in the original skipset.'
    if skipset is None: skipset = set()
    for item in iterable:
        if item not in skipset:
            skipset.add(item)
            yield item


def remove_redundant(metaclasses):
    skipset = set([types.ClassType])
    for meta in metaclasses: # determines the metaclasses to be skipped
        skipset.update(inspect.getmro(meta)[1:])
    return tuple(skip_redundant(metaclasses, skipset))

######### now the core of the module: two mutually recursive functions #########

memorized_metaclasses_map = {}

def get_noconflict_metaclass(bases, left_metas, right_metas):
     '''Not intended to be used outside of this module, unless you know
     what you are doing.'''
     # make tuple of needed metaclasses in specified priority order
     metas = left_metas + tuple(map(type, bases)) + right_metas
     needed_metas = remove_redundant(metas)
     # return existing confict-solving meta, if any
     if needed_metas in memorized_metaclasses_map:
       return memorized_metaclasses_map[needed_metas]
     # nope: compute, memoize and return needed conflict-solving meta
     elif not needed_metas: # wee, a trivial case, happy us
         meta = type
     elif len(needed_metas) == 1: # another trivial case
        meta = needed_metas[0]
     # check for recursion, can happen i.e. for Zope ExtensionClasses
     elif needed_metas == bases:
         raise TypeError('Incompatible root metatypes', needed_metas)
     else: # gotta work ...
         metaname = '_' + ''.join([m.__name__ for m in needed_metas])
         meta = classmaker()(metaname, needed_metas, {})
     memorized_metaclasses_map[needed_metas] = meta
     return meta

def classmaker(left_metas=(), right_metas=()):
    '''
    Avoid conflicts when using metaclass.
    To use this:

    class MyClass(AMetaclassedA, AMetaclassedB):
        __metaclass__ = classmaker()

        # then put the rest of your code ...

    '''
    def make_class(name, bases, adict):
        metaclass = get_noconflict_metaclass(bases, left_metas, right_metas)
        return metaclass(name, bases, adict)
    return make_class

################################################################################




