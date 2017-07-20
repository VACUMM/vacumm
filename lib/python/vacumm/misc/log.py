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
__doc__ = 'Logging features'

import fnmatch, logging, logging.handlers, os, sys, traceback, warnings

try:
    # location since Python 2.7
    from weakref import WeakSet
except ImportError:
    # separately installed
    #warnings.warn('Cannot import weakref.WeakSet, trying weakrefset.WeakSet')
    from weakrefset import WeakSet
    # in debugging mode, show which weakrefset module is in use
    if '--debug' in sys.argv:
        import weakrefset
        warnings.warn('Using %s'%(weakrefset.__file__))

try:
    import curses
    curses.setupterm()
except Exception: hascurses = False
else: hascurses = True

logging.VERBOSE = logging.DEBUG + (logging.INFO - logging.DEBUG) / 2
logging.NOTICE = logging.INFO + (logging.WARNING - logging.INFO) / 2

logging.addLevelName(logging.VERBOSE, 'VERBOSE')
logging.addLevelName(logging.NOTICE, 'NOTICE')

# Patch logging handleError annoying behavior (<2.7?)
# NOTE: logging errors are not just printed, exception is raised
#      (logging.raiseExceptions is used either to mean raise or print ...)
_org_handle_error = logging.Handler.handleError
def handle_error(self, *args, **kwargs):
    if logging.raiseExceptions:
        raise
logging.Handler.handleError = handle_error
del handle_error

class LogLevelError(Exception):
    def __init__(self, level):
        Exception.__init__(self, 'Invalid log level: "%s", available are: %s'%(level, ', '.join(get_str_levels())))

def check_level(level):
    if isinstance(level, basestring):
        level = level.strip().upper()
    if level not in logging._levelNames.keys():
        raise LogLevelError(level)
    return level

def convert_level(level):
    '''
    Convert a level from int/str to str/int
    '''
    if isinstance(level, basestring):
        level = level.strip().upper()
    return logging._levelNames[check_level(level)]

def get_int_levels():
    return sorted(filter(lambda l: isinstance(l, int), logging._levelNames.keys()))

def get_str_levels():
    return map(lambda l: get_str_level(l), get_int_levels())

def get_int_level(level):
    '''
    Get the int level for the given int/str level
    '''
    if isinstance(level, basestring):
        return convert_level(level)
    elif isinstance(level, int):
        return check_level(level)
    else:
        raise LogLevelError(level)

def get_str_level(level):
    '''
    Get the str level for the given int/str level
    '''
    if isinstance(level, basestring):
        return check_level(level)
    elif isinstance(level, int):
        return convert_level(level)
    else:
        raise LogLevelError(level)

def filter_name(names, value, default):
    if not isinstance(names, (list,tuple)):
        names = (names,)
    names = map(lambda n: unicode(n).lower().strip(), names)
    svalue = unicode(value)
    if '=' in svalue:
        for v in svalue.split(','):
            n, v = v.split('=', 1)
            for name in names:
                if fnmatch.fnmatch(name, n.lower().strip()):
                    return v
        return default
    return value


class ColoredStreamHandler(logging.StreamHandler):
    colors = {
        'blackf':30,
        'redf':31,
        'greenf':32,
        'yellowf':33,
        'bluef':34,
        'purplef':35,
        'cyanf':36,
        'whitef':37,
        'blackb':40,
        'redb':41,
        'greenb':42,
        'yellowb':43,
        'blueb':44,
        'purpleb':45,
        'cyanb':46,
        'whiteb':47,
        'boldon':1,
        'boldoff':22,
        'italicson':3,
        'italicsoff':23,
        'ulon':4,
        'uloff':24,
        'invon':7,
        'invoff':27,
        'reset':0,
    }

    level_colors = {
        logging.NOTSET:'reset',
        logging.DEBUG:'cyanf',
        logging.VERBOSE:'bluef',
        logging.INFO:'greenf',
        logging.NOTICE:'greenf,boldon',
        logging.WARNING:'yellowf,boldon',
        logging.ERROR:'redf,boldon',
        logging.CRITICAL:'purplef,boldon',
    }

    def __init__(self, stream, *args, **kwargs):
        '''
        :Params:
            - **kwargs**: may contain:
                - **colorize**: True/False: whether to colorize message if possible (curses)
        '''
        # terminal capabilities (man 5 terminfo)
        # setf: set foreground color
        # setaf: set foreground color using ansi escape
        # setb/setab: same for background
        self.colorable = hascurses and stream.isatty() and (curses.tigetstr('setf') or curses.tigetstr('setaf'))
        self.colorize = kwargs.pop('colorize', True) and self.colorable
        # TODO: add per log record components colorization
        logging.StreamHandler.__init__(self, stream, *args, **kwargs)

    def get_level_color(self, level):
        return self.level_colors[get_int_level(level)]

    def format(self, record):
        formatted = logging.StreamHandler.format(self, record)
        # Check if colorize cannot be done
        if not (self.colorable and self.colorize):
            return formatted
        color = self.get_level_color(record.levelno) + getattr(record, 'color', '')
        if color:
            colors = str(color).strip().lower().split(',')
            formatted = '%s%s\033[%dm'%(
                '\033[%sm'%(';'.join('%s'%self.colors.get(c, self.colors['reset']) for c in colors if c)),
                formatted,
                self.colors['reset'])
        return formatted

class Formatter(logging.Formatter): pass

class _Redirector_(object):
    def __init__(self, func, prefix=''):
        self.func = func
        self.prefix = prefix
    def write(self, buf):
        self.func(self.prefix+buf.rstrip().strip('\n'))
    def flush(self):
        pass

class Logger(logging.getLoggerClass()):
    '''Logger based on the standard python logging package.'''

    initial_level = logging.INFO # The initial default level, which would not be changed by set_default_level
    default_level = initial_level
    default_format = dict(
        default='[%(asctime)s %(name)s %(levelname)s] %(message)s',
        console='[%(name)s %(levelname)s] %(message)s',
        logfile='[%(asctime)s %(name)s %(levelname)s] %(message)s',
    )
    # '%Y-%m-%d %H:%M:%S %Z'
    # '%a %d %b %Y %H:%M:%S %Z'
    default_date_format = dict(
        default='%Y-%m-%d %H:%M:%S %Z',
        console='',#%H:%M:%S',
        logfile='%Y-%m-%d %H:%M:%S %Z',
    )
    default_max_file_size = (2**10) * 500
    default_max_file_backup = 0
    default_encoding = None # 'utf-8'
    default_colorize = True
    shared_handlers = []

    # Store instances of this class, used to configure already created loggers
    # through the apply_class_optparser_options
    # Avoid memory leaks by using WeakSet !
    instances = WeakSet()

    def __init__(self,
            name=None, level=None,
            format=None, date_format=None,
            console=True, colorize=True,
            logfile=None,
            redirect_warnings=False, redirect_stdout=False, redirect_stderr=False,
            name_filters=None,
            config=None):
        '''
        Initialize the logger (level, formats) and its default handlers.

        :Params:
            - **name**: Logger name
            - **level**: Initial logging level name (case insensitive) or number.
            - **format**: logs format
            - **date_format**: logs date format
            - **console**: add a stream handler (see below)
            - **colorize**: colorize terminal stream handlers if possible
            - **logfile**: add a rotating file handler (see below)
            - **redirect_warnings**, optional: Redirect messages issued by :mod:`warnings.warn`.
            - **redirect_stdout**, optional: Redirect messages issued to sys.stdout.
            - **redirect_stderr**, optional: Redirect messages issued to sys.stderr.
            - **name_filters**: for named based configuration, see :meth:`filter_name`
            - **config**: see :meth:`config`

        Available log levels: **debug, verbose, info, notice, warning, error, critical**.

        Console handler, accessible through **self.console_handler**, is enabled
        by default and uses sys.stderr and is created with :meth:`new_rotating_file_handler`.

        Rotating file handler, accessible through **self.file_handler**, is disabled
        by default  and is created with :meth:`new_stream_handler`.
        **Note: for applications with a large number of loggers, you may create the
        file handler with :meth:`new_rotating_file_handler` and pass this handler to
        the constructor of all loggers, this will avoid opening the same file multiple times.**

        The **logfile** and **console** arguments can take the form of:
            - a list or tuple: interpreted as positionnal arguments for the corresponding method.
            - a dict: interpreted as named arguments for the corresponding method.
            - a boolean for stream handler, True to create the handler
            - a string for the file handler, specifying the log file path
            - an existing handler

        For example, you may customize logger creation:

        >>> l = Logger(
                name=os.path.basename(sys.argv[0]),
                level='debug',
                logfile=dict(
                    filePath='output.log', maxFileSize=1024*50, maxFileBkp=1,
                    format='%(asctime)s - %(name)s - %(level)s - %(message)s',
                    date_format='%Y-%m-%d %H:%M:%S %Z',
                    level='info'
                ),
                console=dict(
                    stream=sys.stdout, colorize=False,
                    format='%(name)s - %(levelname)s - %(message)s'
                )
            )

        '''
        if name is None: name = '%s(%s)'%(self.__class__.__name__, id(self))
        if level is None: level = self.get_default_level()
        if format is None: format = self.get_default_format()
        if date_format is None: date_format = self.get_default_date_format()
        if name_filters is None: name_filters = []
        self.__name_filters = name_filters
        self.__logger_class = logging.getLoggerClass()

        # We must set the level both at init and with set_level below
        # because set_level will use get_level if filtered level doesn't match
        self.__logger_class.__init__(self, name, level=self.initial_level)

        self.__format = format.get('default', self.default_format['default']) if isinstance(format, dict) else format
        self.__date_format = date_format.get('default', self.default_date_format['default']) if isinstance(date_format, dict) else date_format
        
        self.skipCaller((self.log, self.verbose, self.notice))

        # Setup stream handler (console)?
        self.console_handler = None
        if console:
            if isinstance(console, logging.Handler):
                self.console_handler = console
            elif isinstance(console, (list, tuple)):
                self.console_handler = self.new_stream_handler(*console)
            else:
                if not isinstance(console, dict):
                    console = {}
                console.setdefault('colorize', colorize if colorize is not None else self.default_colorize)
                console.setdefault('format',  format['console'] if isinstance(format, dict)
                    and 'console' in format else self.default_format['console'])
                console.setdefault('date_format',  date_format['console'] if isinstance(format, dict)
                    and 'console' in date_format else self.default_date_format['console'])
                self.console_handler = self.new_stream_handler(**console)
            self.add_handler(self.console_handler)

        # Setup file handler (logfile) ?
        # This sould be carefully used: 1 file handler = 1 opened file !
        self.file_handler = None
        if logfile:
            if isinstance(logfile, logging.Handler):
                self.file_handler = logfile
            elif isinstance(logfile, (list, tuple)):
                self.file_handler = self.new_rotating_file_handler(*logfile)
            else:
                if isinstance(logfile, basestring):
                    logfile = {'filePath': logfile}
                if not isinstance(logfile, dict):
                    logfile = {}
                logfile.setdefault('format',  format['logfile'] if isinstance(format, dict)
                    and 'logfile' in format else self.default_format['logfile'])
                logfile.setdefault('date_format',  date_format['logfile'] if isinstance(format, dict)
                    and 'logfile' in date_format else self.default_date_format['logfile'])
                self.file_handler = self.new_rotating_file_handler(**logfile)
            self.add_handler(self.file_handler)

        #self.config(name=name, level=level, format=format, date_format=date_format)
        self.set_level(level, filter=True)
        if config: self.config(config)

        # Add global handlers: handlers registered for all Logger
        for hdlr in self.__class__.shared_handlers:
            self.add_handler(hdlr)

        # Store new instance
        self.instances.add(self)

        # Redirections
        self._redirect_warnings = redirect_warnings
        self._redirect_stderr = redirect_stderr
        self._redirect_stdout = redirect_stdout
        self._redirected_warnings = False
        self._redirected_stderr = False
        self._redirected_stdout = False
        self.start_redirections()

    def start_redirections(self, which='all'):
        assert which in ['all', 'stdout', 'stderr', 'warnings']

        if self._redirect_warnings and not self._redirected_warnings:
            self._old_showwarnings = warnings.showwarning
            warnings.showwarning = self.showwarning
            self._redirected_warnings = True

        if self._redirect_stdout and not self._redirected_stdout:
            if not isinstance(self._redirect_stdout, str):
                self._redirect_stdout = 'debug'
            self._old_sys_stdout = sys.stdout
            sys.stdout = _Redirector_(getattr(self, self._redirect_stdout),
                prefix='STDOUT: ')
            self._redirected_stdout = True

        if self._redirect_stderr and not self._redirected_stderr:
            if not isinstance(self._redirect_stderr, str):
                self._redirect_stderr = 'warning'
            self._old_sys_stderr = sys.stderr
            sys.stderr = _Redirector_(getattr(self, self._redirect_stderr),
                prefix='STDERR: ')
            self._redirected_stderr = True

    def stop_redirections(self, which='all'):
        assert which in ['all', 'stdout', 'stderr', 'warnings']

        if (self._redirect_warnings and which in ['all', 'warnings'] and
                 self._redirected_warnings):
            warnings.showwarning = self._old_showwarnings
            self._redirected_warnings = False

        if (self._redirect_stdout and which in ['all', 'stdout'] and
                 self._redirected_stdout):
            sys.stdout = self._old_sys_stdout
            self._redirected_stdout = False
        if (self._redirect_stderr and which in ['all', 'stderr'] and
                 self._redirected_stderr):
            sys.stderr = self._old_sys_stderr
            self._redirected_stderr = False

    def set_name(self, name):
        self.name = name

    def get_name(self):
        return self.name

    def set_name_filters(self, filters):
        self.__name_filters = filters

    def get_name_filters(self):
        return list(self.filters)

    def filter_name(self, value, default):
        '''
        Used by the set_{level,format,date_format} instance methods.
        The logger get_name() and instance variable self.__name_filters are evaluated against the
        couples name1=value1[,name2=value2,...] in value to determine the value to return, default
        is returned if nothing match, value is return if it isn't a filter expression.
        '''
        return filter_name([self.get_name()]+self.__name_filters, value, default)

    def set_level(self, level, filter=False):
        if isinstance(level, dict):
            if 'console' in level and self.console_handler:
                self.console_handler.setLevel(get_int_level(level['console']))
            if 'logfile' in level and self.file_handler:
                self.file_handler.setLevel(get_int_level(level['logfile']))
        else:
            if filter: # warn: '=' character in fmt are not allowed ...
                level = self.filter_name(level, self.get_level_name())
            self.__logger_class.setLevel(self, get_int_level(level))

    setLevel = set_level # do not delete, redefinition of logging.Logger.setLevel

    def set_level_console(self, level, filter=False):
        self.set_level({'console':level}, filter=filter)

    def set_level_logfile(self, level, filter=False):
        self.set_level({'logfile':level}, filter=filter)

    def get_level(self):
        return get_int_level(self.level)

    def get_level_name(self):
        return get_str_level(self.level)

    def set_format(self, format, filter=False):
        if isinstance(format, dict):
            if 'console' in format and self.console_handler:
                self.set_handler_formats(self.console_handler, format=format['console'])
            if 'logfile' in format and self.file_handler:
                self.set_handler_formats(self.file_handler, format=format['logfile'])
        else:
            if filter: # warn: '=' character in fmt are not allowed ...
                format = self.filter_name(format, self.get_format())
            self.__format = format
            for h in self.handlers:
                self.set_handler_formats(h, format=format)

    def get_format(self):
        return self.__format

    def set_date_format(self, format, filter=False):
        if isinstance(format, dict):
            if 'console' in format and self.console_handler:
                self.set_handler_formats(self.console_handler, date_format=format['console'])
            if 'logfile' in format and self.file_handler:
                self.set_handler_formats(self.file_handler, date_format=format['logfile'])
        else:
            if filter: # warn: '=' character in fmt are not allowed ...
                format = self.filter_name(format, self.get_date_format())
            self.__date_format = format
            for h in self.handlers:
                self.set_handler_formats(h, date_format=format)

    def get_date_format(self):
        return self.__date_format

    @classmethod
    def new_formatter(cls, format=None, date_format=None, formatter=None):
        if format is None:
            if formatter and formatter.format:
                format = formatter._fmt
            else:
                format = cls.default_format['default']
        if date_format is None:
            if formatter and formatter.datefmt:
                date_format = formatter.datefmt
            else:
                date_format = cls.default_date_format['default']
        return Formatter(fmt=format, datefmt=date_format)

    @classmethod
    def set_handler_formats(cls, handler, format=None, date_format=None):
        f = cls.new_formatter(format=format, date_format=date_format, formatter=handler.formatter)
        handler.setFormatter(f)
        return f

    def is_verbose(self):
        return self.get_level() <= logging.VERBOSE

    def is_debug(self):
        return self.get_level() <= logging.DEBUG

    def add_handler(self, handler, **kwargs):
        if not handler.formatter:
            handler.setFormatter(self.new_formatter())
        if kwargs.pop('shared', False):
            self.add_shared_handler(handler)
        self.__logger_class.addHandler(self, handler)
    addHandler = add_handler # do not delete, redefinition of logging.Logger.addHandler

    def remove_handler(self, hdlr):
        try: hdlr.close()
        except: print>>sys.stderr, traceback.format_exc()
        try: self.__logger_class.removeHandler(self, hdlr)
        except: print>>sys.stderr, traceback.format_exc()
    removeHandler = remove_handler # do not delete, redefinition of logging.Logger.removeHandler

    def clean(self):
        for hdlr in list(self.handlers):
            self.remove_handler(hdlr)

    @classmethod
    def new_stream_handler(cls, stream=sys.stderr, colorize=False, level=None, format=None, date_format=None, **kwargs):
        handler = ColoredStreamHandler(stream, colorize=colorize)
        if level: handler.setLevel(get_int_level(level))
        if format or date_format: cls.set_handler_formats(handler, format=format, date_format=date_format)
        return handler

    def add_stream_handler(self, *args, **kwargs):
        handler = self.new_stream_handler(*args, **kwargs)
        self.add_handler(handler, **kwargs)
        return handler

    def colorize(self, active=None):
        if self.console_handler and isisntance(self.console_handler, ColoredStreamHandler):
            if active is not None: self.console_handler.colorize = bool(active)
            return self.console_handler.colorize

    @classmethod
    def new_rotating_file_handler(cls, filePath=None, mode='a', maxFileSize=None, maxFileBkp=None, encoding=None, level=None, format=None, date_format=None, **kwargs):
        if maxFileSize is None: maxFileSize = cls.default_max_file_size
        if maxFileBkp is None: maxFileBkp = cls.default_max_file_backup
        if encoding is None: encoding = cls.default_encoding
        if not filePath:
            if sys.argv[0] in ['', '-', '-c']: filePath = 'python'
            else: filePath = os.path.splitext(os.path.basename(sys.argv[0]))[0]
            filePath += '.log'
        fileDir = os.path.dirname(filePath)
        if fileDir and not os.path.isdir(fileDir):
            os.makedirs(fileDir)
        handler = logging.handlers.RotatingFileHandler(filename=filePath, mode=mode, maxBytes=maxFileSize, backupCount=maxFileBkp, encoding=encoding)
        if level: handler.setLevel(get_int_level(level))
        if format or date_format: cls.set_handler_formats(handler, format=format, date_format=date_format)
        return handler

    def add_rotating_file_handler(self, *args, **kwargs):
        handler = self.new_rotating_file_handler(*args, **kwargs)
        self.add_handler(handler, **kwargs)
        return handler

    def notice(self, msg, *args, **kwargs):
        self.log(logging.NOTICE, msg, *args, **kwargs)

    def verbose(self, msg, *args, **kwargs):
        self.log(logging.VERBOSE, msg, *args, **kwargs)

#   def _log(self, level, msg, *args, **kwargs):
#       msg = msg.encode('utf-8', 'replace')
#       return super(self.__class__, self)._log(level, msg, *args, **kwargs)

    def log(self, level, msg, *args, **kwargs):
        level = get_int_level(level)
        return self.__logger_class.log(self, level, msg, *args, **kwargs)

    # redefine findCaller to ignore certain callers,
    # used to support the new levels (verbose/notice)
    def findCaller(self):
        """
        Find the stack frame of the caller so that we can note the source
        file name, line number and function name.
        """
        f = logging.currentframe()
        #On some versions of IronPython, currentframe() returns None if
        #IronPython isn't run with -X:Frames.
        if f is not None:
            f = f.f_back
        rv = "(unknown file)", 0, "(unknown function)"
        skipCallers = getattr(self, '_skipCallers', set())
        while hasattr(f, "f_code"):
            co = f.f_code
            filename = os.path.normcase(co.co_filename)
            skip = False
            if f.f_code in skipCallers:
                skip = True
            if filename == logging._srcfile:
                skip = True
            if skip:
                f = f.f_back
                continue
            rv = (co.co_filename, f.f_lineno, co.co_name)
            break
        return rv
     
    def skipCaller(self, caller=None, skip=True):
        '''
        Register a function, method, code or set/list/tuple of codes to be skipped by findCaller
        '''
        skipCallers = getattr(self, '_skipCallers', set())
        if caller is None:
            return skipCallers
        elif isinstance(caller, (set, list, tuple)):
            for c in caller:
                self.skipCaller(c, skip)
            return
        import types
        if isinstance(caller, types.CodeType):
            func_code = caller
        elif isinstance(caller, types.FunctionType):
            func_code = caller.func_code
        elif isinstance(caller, types.MethodType):
            func_code = caller.__func__.func_code
        else:
            raise TypeError('skipCaller accept only function or method types, found %s'%type(caller))
        if skip:
            skipCallers.add(func_code)
        else:
            skipCallers.remove(func_code)
        setattr(self, '_skipCallers', skipCallers)

    @classmethod
    def set_default_level(cls, level, filter=False):
        if filter: level = filter_name(cls.__name__, level, cls.get_default_level())
        cls.default_level = get_int_level(level)
    @classmethod
    def get_default_level(cls):
        return cls.default_level

    @classmethod
    def set_default_format(cls, format, target=None, filter=False):
        '''
        Set the default format for the specified target.
        The target can be one or a list of: %s.
        If target is None, all default formats are changed.
        '''%(', '.join(cls.default_format))
        if filter: format = filter_name(cls.__name__, format, cls.get_default_format(target)) # warn: '=' character in fmt are then not allowed ...
        if target is None: target = cls.default_format.keys()
        if not isinstance(target, (list, tuple)): target = (target,)
        for k in target: cls.default_format[k] = format

    @classmethod
    def get_default_format(cls, target=None):
        '''Get the default format for the specified target, or all formats dict if target is None'''
        return cls.default_format.copy() if target is None else cls.default_format[target]

    @classmethod
    def set_default_date_format(cls, format, target=None, filter=False):
        '''
        Set the default date format for the specified target.
        The target can be one or a list of: %s.
        If target is None, all default date formats are changed.
        '''%(', '.join(cls.default_date_format))
        if filter: format = filter_name(cls.__name__, format, cls.get_default_date_format(target)) # warn: '=' character in fmt are then not allowed ...
        if target is None: target = cls.default_date_format.keys()
        if not isinstance(target, (list, tuple)): target = (target,)
        for k in target: cls.default_date_format[k] = format

    @classmethod
    def get_default_date_format(cls, target=None):
        '''Get the default date format for the specified target, or all formats dict if target is None'''
        return cls.default_date_format.copy() if target is None else cls.default_date_format[target]

    @classmethod
    def set_default_max_file_size(cls, value, filter=False):
        if filter: value = filter_name(cls.__name__, value, cls.get_default_max_file_size) # warn: '=' character in fmt are then not allowed ...
        cls.default_max_file_size = value

    @classmethod
    def set_default_max_file_backup(cls, value, filter=False):
        if filter: value = filter_name(cls.__name__, value, cls.get_default_max_file_backup) # warn: '=' character in fmt are then not allowed ...
        cls.default_max_file_backup = value

    @classmethod
    def default_colorize(cls, active=None):
        if active is not None: cls.default_colorize = active
        return cls.default_colorize
    get_default_colorize = default_colorize
    set_default_colorize = default_colorize

    @classmethod
    def add_shared_handler(cls, handler):
        cls.shared_handlers.append(handler)

    @classmethod
    def _add_parser_options(cls, parser, rename=None, prefix='', suffix=''):
        '''Add logging options to a optparse OptionParser'''
        names = dict(((n,n) for n in ('debug', 'verbose', 'loglevel','logformat','logdateformat','lognocolor')))
        if rename is not None: names.update(rename)
        fmt = lambda n: '--%s%s%s'%(prefix, names[n], suffix)
        import optparse
        # we must handle both optparse and argparse
        add = getattr(parser, 'add_option' if isinstance(parser, optparse.OptionParser) else 'add_argument')
        add(fmt('debug'), dest='_logger_debug', action='store_true', default=False, help='Set logging level to DEBUG')
        add(fmt('verbose'), dest='_logger_verbose', action='store_true', default=False, help='Set logging level to VERBOSE')
        add(fmt('loglevel'), dest='_logger_level', action='store', default=None, metavar='level', help='Set logging level, available levels are: %s. You may restrict this level for specific classes with filters: --loglevel=Class1=debug,Class2=warning'%(', '.join(get_str_levels())))
        add(fmt('logformat'), dest='_logger_format', action='store', default=None, metavar='format', help='Set logging format. Can also use filters.')
        add(fmt('logdateformat'), dest='_logger_dateformat', action='store', default=None, metavar='format', help='Set logging date format. Can also use filters.')
        add(fmt('lognocolor'), dest='_logger_nocolor', action='store_true', default=None, help='Disable logging colors')

    @classmethod
    def add_optparser_options(cls, *args, **kwargs):
        return cls._add_parser_options(*args, **kwargs)

    @classmethod
    def add_argparser_options(cls, *args, **kwargs):
        return cls._add_parser_options(*args, **kwargs)

    @classmethod
    def _apply_class_parser_options(cls, options, instances=True):
        '''Apply logging options from a optparse OptionParser, to this class defaults and optionally all existing instances'''
        if options._logger_debug: cls.set_default_level(logging.DEBUG)
        if options._logger_verbose: cls.set_default_level(logging.VERBOSE)
        if options._logger_level: cls.set_default_level(options._logger_level)
        if options._logger_format: cls.set_default_format(options._logger_format)
        if options._logger_dateformat: cls.set_default_date_format(options._logger_dateformat)
        if options._logger_nocolor: cls.set_default_colorize(False)
        if instances:
            for logger in cls.instances:
                logger.apply_optparser_options(options)

    @classmethod
    def apply_class_optparser_options(cls, *args, **kwargs):
        return cls._apply_class_parser_options(*args, **kwargs)

    @classmethod
    def apply_class_argparser_options(cls, *args, **kwargs):
        return cls._apply_class_parser_options(*args, **kwargs)

    def _apply_parser_options(self, options):
        '''Apply logging options from a optparse OptionParser, to this instance'''
        if options._logger_debug: self.set_level(logging.DEBUG)
        if options._logger_verbose: self.set_level(logging.VERBOSE)
        if options._logger_level: self.set_level(options._logger_level, filter=True)
        if options._logger_format: self.set_format(options._logger_format, filter=True)
        if options._logger_dateformat: self.set_date_format(options._logger_dateformat, filter=True)
        if options._logger_nocolor: self.colorize(False)

    def apply_optparser_options(self, *args, **kwargs):
        return self._apply_parser_options(*args, **kwargs)

    def apply_argparser_options(self, *args, **kwargs):
        return self._apply_parser_options(*args, **kwargs)

    def config(self, *args, **kwargs):
        '''
        Get or set logging configuration

        This method uses the get_* and set_* methods.

        :Params:
            - **args**: depending of the type:
                - basestring: (getter) return the corresponding value(s)
                - Logger: (setter) use this logger's configuration
                - dict: (setter) use this dict's configuration
            - **kwargs**: key/value pairs of configuration which are to be set

        Getter:

            >>> logger.config('level')
            %(default_level)r
            >>> logger.config('level', 'format')
            (%(default_level)r, %(default_format)r)
            >>> logger.config() # return complete configuration as a dict
            TODO

        Setter:

            >>> logger.config(level=%(default_level)r)
            >>> logger.config(level=%(default_level)r, format=%(default_format)r)
            >>> logger.config(configdict) # use a dict as config key/value pairs
            >>> logger.config(otherlogger) # use another logger's configuration

        '''
        args = list(args)
        if kwargs:
            args.insert(0, kwargs)
        if not args:
            # 'name' should not appear in this list
            args = ['level', 'format', 'date_format']
        c = {}
        for arg in args:
            if isinstance(arg, basestring):
                c[arg] = getattr(self, 'get_'+arg)()
            elif isinstance(arg, dict):
                for k,v in arg.iteritems():
                    getattr(self, 'set_'+k)(v)
            elif isinstance(arg, Logger):
                self.config(arg.config())
            else:
                raise ValueError('Cannot configure using objects of type %r'%(arg))
        return c

    # TODO: merge with config method above
    def load_config(self, cfgo=None, nested=True):
        if cfgo is None: return
        if nested:
            cfgo = cfgo.get(isinstance(nested, basestring) and nested or self.__class__.__name__, cfgo)
        if cfgo.get('level', None): self.set_level(cfgo['level'])
        if cfgo.get('format', None): self.set_format(cfgo['format'])
        if cfgo.get('date_format', None): self.set_date_format(cfgo['date_format'])
        if cfgo.get('file', None): self.add_rotating_file_handler(filePath=cfgo['file'])

    def showwarning(self, message, category, filename, lineno,
            file=None):
        self.warning(
            'REDIRECTED: %s:%s: %s:%s',
            filename, lineno,
            category.__name__, message,
        )


#   @classmethod
#   def load_default_config(cls, *args, **kwargs):
#       cls.load_config(cls, *args, **kwargs)


#class ClassFilter(logging.Filter):
#   def __init__(self, obj):
#       self.obj = obj
#   def add_rule(self, rule):
#       self.rule.append(rule)
#       self.eval_rule()
#   def eval_rule():
#       self.filtered = False
#   def filter(self, record):
#       return self.filtered


# Fix logger adapter is not a new-style class in python < 2.7, resulting
# in super() failures when inheriting LoggerAdapter
if sys.version_info[:2] < (2,7):
    class _LoggerAdapter(object, logging.LoggerAdapter):
        def __init__(self, *a, **k):
            object.__init__(self)
            logging.LoggerAdapter.__init__(self, *a, **k)
else:
   _LoggerAdapter = logging.LoggerAdapter

class LoggerAdapter(_LoggerAdapter):
    def __init__(self, **kwargs):
        kwargs.setdefault('name', self.__class__.__name__)
        logger = Logger(**kwargs)
        _LoggerAdapter.__init__(self, logger, {})

    def get_logger(self):
        return self.logger

    def process(self, msg, kwargs):
        import copy
        kwargs = copy.deepcopy(kwargs)
        extra = kwargs.pop('extra', {})
        msg, kwargs = logging.LoggerAdapter.process(self, msg, kwargs)
        kwargs = copy.deepcopy(kwargs)
        kwargs['extra'].update(extra)
        return msg, kwargs

    def verbose(self, msg, *args, **kwargs):
        """
        Delegate a log call to the underlying logger, after adding
        contextual information from this adapter instance.
        """
        msg, kwargs = self.process(msg, kwargs)
        self.logger.verbose(msg, *args, **kwargs)

    def notice(self, msg, *args, **kwargs):
        """
        Delegate a log call to the underlying logger, after adding
        contextual information from this adapter instance.
        """
        msg, kwargs = self.process(msg, kwargs)
        self.logger.notice(msg, *args, **kwargs)

    def is_debug(self):
        '''Return True if the logger level is debug (or greater)'''
        return self.logger.is_debug()

    def is_verbose(self):
        '''Return True if the logger level is verbose (or greater)'''
        return self.logger.is_verbose()

#
# MAY HAVE SIDE EFFECTS WHEN DAEMONIZED, TOCHECK: DO NOT CLOSE STREAMS WHEN DAEMONIZING
# bug can be reproduced when e.g. MySQLdb send warnings
# Redirect stdout/stderr to a logger
class StreamLogWrapper(object):
    def __init__(self, stream, logger, name=''):
        self.stream = stream
        self.logger = logger
        self.name = name
    def __getattr__(self, a):
        return getattr(self.stream, a)
    def write(self, s):
        if s.endswith('\n'): s = s[:-1]
        if not s: return
        self.logger.warning('%s: %s', self.name, s)

# A default logger
default = logger = Logger.default = Logger(name=os.path.basename(sys.argv[0] in ('', '-c') and 'python' or sys.argv[0]))
if '--verbose' in sys.argv: default.set_level('verbose')
if '--debug' in sys.argv: default.set_level('debug')
# Register logging function using the default logger
for l in map(str.lower, get_str_levels()):
    f = getattr(default, l, None)
    if callable(f): globals()[l] = f
del l
# Plus the special exception log func
exception = default.exception


def _testLog():

    levels = []
    for i in sorted(filter(lambda l: isinstance(l, int), logging._levelNames.keys())):
        levels.append(i)
        levels.append(get_str_level(i))
    funcs = filter(lambda l: isinstance(l, basestring) and l != 'NOTSET', levels)
    for lvl in levels:
        log = Logger(name='Logger(level=%s)'%lvl, level=lvl, format='%(name)-30s : message sent with level %(levelname)-10s: %(message)s')
        for f in funcs:
            getattr(log, f.lower())('Hello')

#   Logger.parseOptions()
    class A(Logger): pass
    a = A()
    for f in funcs:
        getattr(a, f.lower())('Hello log function %s'%f)

    class B(Logger):
        def __init__(self, **kwargs):
            pfx = 'log_'
            lkw = dict(map(lambda (k,v): (k[len(pfx):],v) , filter(lambda (k,v): k.startswith(pfx), kwargs.iteritems())))
            Logger.__init__(self, **lkw)
    B(log_level='notset').debug('Hello')
    B(log_level='info').debug('Hello')
    B(log_format='%(message)s').info('Hello')


if __name__ == '__main__':
    _testLog()


