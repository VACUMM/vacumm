#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar (2010)
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
    warnings.warn('Cannot import weakref.WeakSet, trying weakrefset.WeakSet')
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
#       (logging.raiseExceptions is used either to mean raise or print ...)
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
        # terminal capabilities (man 5 terminfo)
        # setf: set foreground color
        # setaf: set foreground color using ansi escape
        # setb/setab: same for background
        self.colorable = hascurses and stream.isatty() and (curses.tigetstr('setf') or curses.tigetstr('setaf'))
        #####
        if os.environ.get('LOG_FORCE_COLORIZE', None): self.colorable = True
        #####
        self.colorize = self.colorable and True
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
                #''.join('\033[%dm'%self.colors.get(c, self.colors['reset']) for c in colors if c),
                '\033[%sm'%(';'.join('%s'%self.colors.get(c, self.colors['reset']) for c in colors if c)),
                formatted,
                self.colors['reset'])
        return formatted

class Formatter(logging.Formatter): pass

class Logger(logging.getLoggerClass()):
    '''Logger based on the standard python logging package.
    
    :Params:
        - **name**: Logger name
        - **level**: Initial logging level, case insensitive. One of:
            - debug
            - verbose
            - info
            - notice
            - warning
            - error
            - critical
        - **format**: logs format
        - **date_format**: logs date format
        - **console**: add a stream handler using sys.stderr
        - **colorize**: colorize terminal stream handlers if possible
        - **name_filters**: optional, for named based configuration, see :meth:`filter_name`
        
    '''
    initial_level = logging.INFO # The initial default level, which would not be changed by set_default_level
    default_level = initial_level
    default_format = '[%(asctime)s %(name)s %(levelname)s] %(message)s'
    default_date_format = '%Y-%m-%d %H:%M:%S %Z' # '%a %d %b %Y %H:%M:%S %Z'
    default_max_file_size = (2**10) * 500
    default_max_file_backup = 0
    default_encoding = None # 'utf-8'
    default_colorize = True
    shared_handlers = []
    
    # Avoid memory leaks !
    instances = WeakSet()
    
    def __init__(self, name=None, level=None,
            format=None, date_format=None,
            console=True, colorize=None,
            name_filters=None):
        if name is None: name = '%s(%s)'%(self.__class__.__name__, id(self))
        if level is None: level = self.default_level
        if format is None: format = self.default_format
        if date_format is None: date_format = self.default_date_format
        if colorize is None: colorize = self.default_colorize
        if name_filters is None: name_filters = []
        self.name_filters = name_filters
        self.__logger_class = logging.getLoggerClass()
        # We must set the level both at init and with set_level below
        # because set_level will use get_level if filtered level doesn't match
        self.__logger_class.__init__(self, name, level=self.initial_level)
        self.__format = format
        self.__date_format = date_format
        self.__formatter = Formatter(fmt=self.__format, datefmt=self.__date_format)
        self.set_level(level, filter=True)
        # If stderr is True, define a default handler to sys.stderr
        self.console_handler = None
        if console: self.console_handler = self.add_stream_handler(colorize=colorize)
        # Add global handlers: handlers registered for all Logger
        for hdlr in self.__class__.shared_handlers:
            self.add_handler(hdlr)
        self.instances.add(self)
    
    # PARTIAL IMPLEMENTATION: handlers also have a level check ...!
#    def log(self, level, msg, *args, **kwargs):
#        ''' Add a keyword argument 'force' to bypass the isEnabledFor(level) check.
#        '''
#        if not isinstance(level, int):
#            if logging.raiseExceptions:
#                raise TypeError('level must be an integer')
#            else:
#                return
#        if self.isEnabledFor(level) or kwargs.get('force', False):
#            self._log(level, msg, args, **kwargs)
    
    def set_name(self, name): self.name = name
    def get_name(self): return self.name
    
    def filter_name(self, value, default):
        '''
        Used by the set_{level,format,date_format} instance methods.
        The logger get_name() and instance variable self.name_filters are evaluated against the
        couples name1=value1[,name2=value2,...] in value to determine the value to return, default
        is returned if nothing match, value is return if it isn't a filter expression.
        '''
        return filter_name([self.get_name()]+self.name_filters, value, default)
    
    def set_level(self, level, filter=False):
        if filter: level = self.filter_name(level, self.get_level_name())
        self.__logger_class.setLevel(self, get_int_level(level))
    setLevel = set_level # do not delete, redefinition of logging.Logger.setLevel
    def get_level(self): return get_int_level(self.level)
    def get_level_name(self): return get_str_level(self.level)
    
    def set_format(self, format, filter=False):
        if filter: format = self.filter_name(format, self.get_format())
        self.__format = format
        self.__formatter = Formatter(fmt=self.__format, datefmt=self.__date_format)
        for h in self.handlers:
            h.setFormatter(self.__formatter)
    def get_format(self): return self.__format
    
    def set_date_format(self, date_format, filter=False):
        # warn: '=' character in fmt are then not allowed ...
        if filter: format = self.filter_name(format, self.get_date_format())
        self.__date_format = date_format
        self.__formatter = Formatter(fmt=self.__format, datefmt=self.__date_format)
        for h in self.handlers:
            h.setFormatter(self.__formatter)
    def get_date_format(self): return self.__date_format
    
    def is_verbose(self):
        return self.get_level() <= logging.VERBOSE
    def is_debug(self):
        return self.get_level() <= logging.DEBUG
    
    def add_handler(self, handler, **kwargs):
        if not handler.formatter:
            handler.setFormatter(self.__formatter)
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
    def new_stream_handler(cls, stream=sys.stderr, colorize=False, **kwargs):
        return (ColoredStreamHandler if colorize else logging.StreamHandler)(stream)
    
    def add_stream_handler(self, *args, **kwargs):
        handler = self.new_stream_handler(*args, **kwargs)
        self.add_handler(handler, **kwargs)
        return handler
    
    def colorize(self, active=None):
        if self.console_handler:
            if active is not None: self.console_handler.colorize = bool(active)
            return self.console_handler.colorize
        
    @classmethod
    def new_rotating_file_handler(cls, filePath=None, mode='a', maxFileSize=None, maxFileBkp=None, encoding=None, **kwargs):
        if maxFileSize is None: maxFileSize = cls.default_max_file_size
        if maxFileBkp is None: maxFileBkp = cls.default_max_file_backup
        if encoding is None: encoding = cls.default_encoding
        if not filePath:
            filePath = self.name+'.log'
        fileDir = os.path.dirname(filePath)
        if fileDir and not os.path.isdir(fileDir):
            os.makedirs(fileDir)
        return logging.handlers.RotatingFileHandler(filename=filePath, mode=mode, maxBytes=maxFileSize, backupCount=maxFileBkp, encoding=encoding)
    
    def add_rotating_file_handler(self, *args, **kwargs):
        handler = self.new_rotating_file_handler(*args, **kwargs)
        self.add_handler(handler, **kwargs)
        return handler
    
    def notice(self, msg, *args, **kwargs):
        self.log(logging.NOTICE, msg, *args, **kwargs)
    
    def verbose(self, msg, *args, **kwargs):
        self.log(logging.VERBOSE, msg, *args, **kwargs)
    
#    def _log(self, level, msg, *args, **kwargs):
#        msg = msg.encode('utf-8', 'replace')
#        return super(self.__class__, self)._log(level, msg, *args, **kwargs)
    
    @classmethod
    def set_default_level(cls, level, filter=False):
        if filter: level = filter_name(cls.__name__, level, cls.get_default_level())
        cls.default_level = get_int_level(level)
    @classmethod
    def get_default_level(cls):
        return cls.default_level
    
    @classmethod
    def set_default_format(cls, format, filter=False):
        if filter: format = filter_name(cls.__name__, format, cls.get_default_format()) # warn: '=' character in fmt are then not allowed ...
        cls.default_format = format
    @classmethod
    def get_default_format(cls):
        return cls.default_format
    
    @classmethod
    def set_default_date_format(cls, format, filter=False):
        if filter: format = filter_name(cls.__name__, format, cls.get_default_date_format()) # warn: '=' character in fmt are then not allowed ...
        cls.default_date_format = format
    @classmethod
    def get_default_date_format(cls):
        return cls.default_date_format
    
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
    def add_optparser_options(cls, parser, rename=None, prefix='', suffix=''):
        names = dict(((n,n) for n in ('debug', 'verbose', 'loglevel','logformat','logdateformat','lognocolor')))
        if rename is not None: names.update(rename)
        fmt = lambda n: '--%s%s%s'%(prefix, names[n], suffix)
        parser.add_option(fmt('debug'), dest='_logger_debug', action='store_true', default=False, help='Set logging level to DEBUG')
        parser.add_option(fmt('verbose'), dest='_logger_verbose', action='store_true', default=False, help='Set logging level to VERBOSE')
        parser.add_option(fmt('loglevel'), dest='_logger_level', action='store', default=None, metavar='level', help='Set logging level, available levels are: %s. You may restrict this level for specific classes with filters: --loglevel=Class1=debug,Class2=warning'%(', '.join(get_str_levels())))
        parser.add_option(fmt('logformat'), dest='_logger_format', action='store', default=None, metavar='format', help='Set logging format. Can also use filters.')
        parser.add_option(fmt('logdateformat'), dest='_logger_dateformat', action='store', default=None, metavar='format', help='Set logging date format. Can also use filters.')
        parser.add_option(fmt('lognocolor'), dest='_logger_nocolor', action='store_true', default=None, metavar='format', help='Disable logging colors')
    
    def apply_optparser_options(self, options):
        if options._logger_debug: self.set_level(logging.DEBUG)
        if options._logger_verbose: self.set_level(logging.VERBOSE)
        if options._logger_level: self.set_level(options._logger_level, filter=True)
        if options._logger_format: self.set_format(options._logger_format, filter=True)
        if options._logger_dateformat: self.set_date_format(options._logger_dateformat, filter=True)
        if options._logger_nocolor: self.colorize(False)
    
    @classmethod
    def apply_class_optparser_options(cls, options, instances=True):
        if options._logger_debug: cls.set_default_level(logging.DEBUG)
        if options._logger_verbose: cls.set_default_level(logging.VERBOSE)
        if options._logger_level: cls.set_default_level(options._logger_level)
        if options._logger_format: cls.set_default_format(options._logger_format)
        if options._logger_dateformat: cls.set_default_date_format(options._logger_dateformat)
        if options._logger_nocolor: cls.set_default_colorize(False)
        if instances:
            for logger in cls.instances:
                logger.apply_optparser_options(options)
    
    def load_config(self, cfgo=None, nested=True):
        if cfgo is None: return
        if nested:
            cfgo = cfgo.get(isinstance(nested, basestring) and nested or self.__class__.__name__, cfgo)
        if cfgo.get('level', None): self.set_level(cfgo['level'])
        if cfgo.get('format', None): self.set_format(cfgo['format'])
        if cfgo.get('date_format', None): self.set_date_format(cfgo['date_format'])
        if cfgo.get('file', None): self.add_rotating_file_handler(filePath=cfgo['file'])
    
#    @classmethod
#    def load_default_config(cls, *args, **kwargs):
#        cls.load_config(cls, *args, **kwargs)


#class ClassFilter(logging.Filter):
#    def __init__(self, obj):
#        self.obj = obj
#    def add_rule(self, rule):
#        self.rule.append(rule)
#        self.eval_rule()
#    def eval_rule():
#        self.filtered = False
#    def filter(self, record):
#        return self.filtered


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
    
    # TODO: fix funcName
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
for l in get_str_levels():
    globals()[l.lower()] = getattr(default, l.lower(), None)
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
    
#    Logger.parseOptions()
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


