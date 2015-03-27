# -*- coding: utf-8 -*-

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

from __future__ import absolute_import
if __name__ == '__main__':
    import sys; del sys.path[0]

__author__ = 'Stéphane Raynaud, Jonathan Wilkins'
__email__ = 'raynaud@actimar.fr, wilkins@actimar.fr'
__doc__ = '''\
Configuration management with validation capabilities and commandline interaction.
It is based on :mod:`configobj` and :mod:`validate` modules.
'''

import copy, datetime, inspect, os, operator, re, sys, shutil, traceback
from optparse import OptionParser, OptionGroup, OptionContainer, Option,  OptionValueError, IndentedHelpFormatter, _
from argparse import (ArgumentParser, _ArgumentGroup, HelpFormatter as ArgHelpFormatter,
    _StoreTrueAction as AP_StoreTrueAction,
    _HelpAction, Action as AP_Action, SUPPRESS as AP_SUPPRESS)
from warnings import warn

from configobj import ConfigObj, flatten_errors
from validate import Validator, ValidateError, VdtTypeError, VdtValueTooSmallError, VdtValueTooBigError
import validate

try: import numpy
except: numpy = None

__all__ = ['ConfigException', 'ValidationWarning', 'ConfigManager', 'print_short_help',
    'opt2rst', 'cfg2rst', 'cfgargparse', 'cfgoptparse', 'getspec', 'get_secnames',
    ]

class ConfigException(Exception):
    pass

class ValidationWarning(Warning):
    pass

class VdtSizeError(ValidateError):
    '''List size is incorrect (nmin, nmax, odd, even or shape)'''
    pass

def _valwrap_(validator):
    '''
    Wrap a validation function to allow extraneous named arguments in specfile,
    this is usefull when getting specification with
    validator._parse_with_caching(configspec[section][option])
    '''
    def validator_wrapper(value, *args, **kwargs):
        # Handle None and default
#        if str(value) == 'None': return None
#        default = args[0] if len(args) else kwargs.get('default', None)
#        if value == '': return default
        # Remove extraneous arguments the validator can't handle
        argspec = inspect.getargspec(validator)
        kwargs = kwargs.copy()
        for k in kwargs.keys():
            if k not in argspec.args:
                kwargs.pop(k)
        return validator(value, *args, **kwargs)
    validator_wrapper.__name__ += '-'+validator.__name__
    return validator_wrapper

def _valwraplist_(validator):
    '''
    Wrap a validation function to handle list value using an existing validator.
    This adds the following list checking behaviors that can be used as named
    arguments in specifications:
        - **n**: required fixed number of elements
        - **nmin**: required minimum number of elements
        - **nmax**: required maximum number of elements
        - **odd**: number of elements must be odd
        - **even**: number of elements must be even
        - **shape**: check value shape (requires numpy)
    '''
    def list_validator_wrapper(value, *args, **kwargs):
        # Handle None and default
        if str(value) == 'None': return None
        default = args[0] if len(args) else kwargs.get('default', ())
        if value == '': value = default
        # Handle single value
        #if not isinstance(value, (list, tuple)): value = (value,)
        # Do list checks
        n = kwargs.pop('n', None)
        if n is not None and len(value) != int(n): raise VdtSizeError('Incorrect size: %s, %s values expected'%(len(value), n))
        nmin = kwargs.pop('nmin', None)
        if nmin is not None and len(value) < int(nmin): raise VdtSizeError('Incorrect size: %s, at least %s values expected'%(len(value), n))
        nmax = kwargs.pop('nmax', None)
        if nmax is not None and len(value) > int(nmax): raise VdtSizeError('Incorrect size: %s, at most %s values expected'%(len(value), nmax))
        odd = validate.is_boolean(kwargs.pop('odd', False))
        if odd and not len(value) % 2: raise VdtSizeError('Incorrect size: %s, odd number of values expected'%(len(value)))
        even = validate.is_boolean(kwargs.pop('even', False))
        if even and len(value) % 2: raise VdtSizeError('Incorrect size: %s, even number of values expected'%(len(value)))
        # TODO?: use numpy to also check min and max ? (applicable on number only...)
        shape = kwargs.pop('shape', None)
        if shape is not None:
            if numpy is None:
                warn('Cannot check shape, numpy package is missing')
            else:
                try:
                    shape, vshape = map(int, shape), numpy.shape(value)
                    if vshape != shape:
                        raise VdtSizeError('Incorrect shape: %s, %s shape expected'%(vshape, shape))
                except Exception, e:
                    raise ValidateError('Cannot test value shape, this may be caused by an irregular array-like shape. Error was:\n%s'%traceback.format_exc())
        # Preserve tuple type
        istuple = isinstance(value, tuple)
        # Validate each values
        value = map(lambda v: validator(v, *args[1:], **kwargs), value)
        return tuple(value) if istuple else value
    list_validator_wrapper.__name__ += '-'+validator.__name__
    return list_validator_wrapper

def _validator_minmax_(value, min=None, max=None, default=(0, 100), type='float'):
    """Validator of a min,max pair"""
    if str(value)=='None': return None
    if isinstance(value, list):
        value = tuple(value)
    elif isinstance(value, (str, unicode)):
        tmp=eval(value)
        try:
            value = eval(value)
        except:
            raise VdtTypeError(value)
    if not len(value)==2: raise VdtTypeError(value)
    out = ()
    for val in value:
        if isinstance(val, (str, unicode)):
            try:
                val = float(val)
            except:
                raise VdtTypeError(value)
        if min is not None and val<float(min): val = float(min)
        if max is not None and val>float(max): val = float(max)
        if type=='int':
            val = int(val)
        out += val,
    if out[0]>out[1]: return out[::-1]
    return out

def _validator_bbox_(value, default=None):
    '''Parse bbox coordinates with value format: x1,y1,x2,y2'''
    if str(value) == 'None': return None
    if value == '': value = default
    c = []
    # two possible delimiters: whitespaces and ','
    for v in value.split(): c.extend(v.split(','))
    if len(c) != 4: raise VdtTypeError(value)
    return map(float, c)

def _validator_figsize_(value, default=(6, 6)):
    """Validator of a figure size (xsize,ysize)"""
    if isinstance(value, list):
        value = tuple(value)
    elif isinstance(value, (str, unicode)):
        try:
            value = eval(value)
        except:
            raise VdtTypeError(value)
    if not len(value)==2: raise VdtTypeError(value)
    out = ()
    for val in value:
        if isinstance(val, (str, unicode)):
            try:
                val = float(val)
            except:
                raise VdtTypeError(value)
        val = max(0, val)
        val = min(20, val)
        out += val,
    return out

def _validator_interval_(value, default=None):
    """Validator of an interval of coordinates (min, max [,bounds])"""
    if str(value) == 'None': return None
    if not isinstance(value, basestring):
        if not isinstance(value, list):
            raise VdtTypeError(value)
        value = ','.join(value)
    if value.startswith('('): value = value[1:]
    if value.endswith(')'): value = value[:-1]
    values = value.split(',')
    if len(values)<2 or len(values)>3: raise VdtTypeError(value)
    out = ()
    for val in values[:2]:
        try:
            val = eval(val)
        except:
            pass
        out += val,
    if len(values)==3 and values[2]:
        m =  re.search('([co]{1,2}[ne]{0,2})', values[2])
        if m is None:
            raise VdtTypeError(value)
        out += m.group(1),
    return out

class VdtDateTimeError(ValidateError):
    pass

def _validator_datetime_(value, default=None, fmt='%Y-%m-%dT%H:%M:%S'):
    '''Parse value as a traditionnal python datetime object'''
    if str(value) == 'None': return None
    if value == '': value = default
    try: return datetime.datetime.strptime(value, fmt)
    except ValueError, e: raise VdtDateTimeError(e)

_re_validator_workdir_split_ = re.compile('\s*>\s*').split
_re_validator_workdir_start_ =  re.compile('^[\[(]').search
_re_validator_workdir_stop_ =  re.compile('[\])]').search
def _validator_workdir_(value, default=''):
    if str(value) == 'None': return ''
    value = value.strip()
    start = _re_validator_workdir_start_(value) is not None
    stop = _re_validator_workdir_stop_(value) is not None
    if (start and not stop) or (stop and not start):
        raise VdtTypeError(value)
    svalue = _re_validator_workdir_split_(value.strip('[]'))
    if len(svalue)==3:
        return '(%s>%s)%s'%tuple(svalue)
    if len(svalue)==2:
        return '(%s>%s)'%tuple(svalue)
    return svalue[0]

# TODO: fix interpolation and expand !
def _validator_path_(value, default='', expand=None):
    '''
    Parse a value as a path
    :Params:
        -- **expand**: expandvars and expandhome on loaded path

    **Warning: expand currently can't work with interpolation**
    '''
    if str(value) == 'None': return None
    if value == '': value = default
    if expand and isinstance(value, basestring):
        return os.path.expandvars(os.path.expanduser(value))
    return value

def _validator_cdtime_(value, min=None, max=None, default=None):
    """Validator of a date (compatible with :func:`cdtime.s2c`)"""
    import cdtime
    value = str(value).strip()
    if not value[0].isdigit(): return value.upper()
    try:
        val = cdtime.s2c(value)
    except:
        raise VdtTypeError(value)
    #not re.match('^\d+(-\d+(-\d+( \d+:(\d+(:\d+(\.(\d+)?)?)?)?)?)?)?$', value):
    if min is not None and val<cdtime.s2c(min): raise VdtValueTooSmallError(value)
    if max is not None and val>cdtime.s2c(max): raise VdtValueTooBigError(value)
    return value

# Define additionnal specifications
# Value should be dict for internal use of this module (iterable, opttype, ...)
# If value is not a dict, it is supposed to be the validator function
_validator_specs_ = {
        # copy of some validate.Validator.functions to later build plural forms
        'integer':validate.is_integer,
        'float':validate.is_float,
        'boolean':validate.is_boolean,
        'string':validate.is_string,
        # single value
        'date':_validator_cdtime_,
        'minmax':_validator_minmax_,
        'figsize':_validator_figsize_,
        'workdir':_validator_workdir_,
        'bbox':_validator_bbox_,
        'datetime':_validator_datetime_,
        'file':_validator_path_,
        'directory':_validator_path_,
        'interval':_validator_interval_,
        # lists validators for these scalars will be automatically generated
}

# 1. Fix specs dicts
# 2. Generate list validators
for k, v in _validator_specs_.items():
    # Check type of spec
    if not isinstance(v, dict):
        v = dict(func=v)
        # Update specs mapping
        _validator_specs_[k] = v
    # Check minimum settings
    v.setdefault('func', validate.is_string)
    v['func'] = _valwrap_(v['func'])
    v.setdefault('iterable', False)
    v.setdefault('opttype', k)
    v.setdefault('argtype', v['func'])
    # Add plural forms and validators to handle list values
    if k.endswith('y'): nk = k[:-1]+'ies'
    elif k.endswith('x'): nk = k+'es'
    else: nk = k+'s'
    if nk not in _validator_specs_:
        nv = v.copy()
        nv['func'] = _valwraplist_(v['func'])
        nv['iterable'] = True
        # OptionParser will check each value, not the list thus we provide the value validator
        nv['opttype'] = nk
        _validator_specs_[nk] = nv

# Build the mapping suitable for Validator.functions
_validator_functions_ = dict((k, v['func']) for k,v in _validator_specs_.iteritems() if 'func' in v)

class ConfigManager(object):
    """A configuration management class based on a configuration specification file
    and a :class:`validate.Validator`

    :Example:

        >>> Cfg = Config('config.ini', interpolation='template')
        >>> Cfg.opt_parse()
        >>> cfg = Cfg.load('config.cfg')

    :See also: :class:`configobj.ConfigObj` and :class:`validate.Validator`

    """
    def __init__(self, cfgspecfile=None, validator=None, interpolation='template',
        encoding=None, boolean_false=True, splitsecdesc=False):
        '''
        :Params:
            - **cfgspecfile**, optional: The specification file to be used with this.
            - **validator**, optional: A custom :class:`validate.Validator`
              to use or a mapping dict of validator functions.
            - **interpolation**, optional: See :class:`configobj.ConfigObj`.
            - **boolean_false**, optional: Make sure that booleans have a default value.
            - **splitsecdesc**, optional: Section descriptions are split in two
              components separated by ':'.
        '''
        # Specifications
        self._encoding = encoding
        if cfgspecfile is not None and isinstance(cfgspecfile, basestring) and not os.path.exists(cfgspecfile):
            raise ConfigException('Specification file not found: %s'%cfgspecfile)
        self._configspecfile = cfgspecfile
        if isinstance(cfgspecfile, ConfigObj):
            self._configspec = cfgspecfile
        else:
            self._configspec = ConfigObj(cfgspecfile, list_values=False,
                interpolation=False, encoding=encoding, raise_errors=True)

        # Validator
        if isinstance(validator, Validator):
            self._validator = validator
        else:
            self._validator = Validator()
        # Wrap default functions to handle none and extra args
        for k, v in self._validator.functions.items():
            self._validator.functions[k] = _valwrap_(v)
        # This modules's validator functions are already wrapped
        self._validator.functions.update(_validator_functions_)
        # If the passed validator is a function mapping
        if isinstance(validator, dict):
            # Also wrap this mapping functions
            validator = validator.copy()
            for k, v in validator.items():
                validator[k] = _valwrap_(v)
            self._validator.functions.update(validator)

        # Makes sure that booleans have a default value
        self._boolean_false = boolean_false
        if boolean_false:
            self._configspec.walk(_walker_set_boolean_false_by_default_,
                validator=self._validator)

        # Interpolation
        if interpolation is True: interpolation='template'
        self._interpolation = interpolation

    def get_spec(self, sec, key):
        '''
        See :func:`getspec`

        If sec is a basestring, use configspec[sec][key]
        Otherwise use sec as a configspec, sec[key]
        '''
        return getspec((self._configspec[sec] if isinstance(sec, basestring) else sec)[key], validator=self._validator)
    getspec = get_spec

    def defaults(self, nocomments=False, interpolation=None):
        """Get the default config

        :Params:

            - **nocomments**, optional: Do not include option comments in config file.
              If equal to 2, remove section comments too.
            - **interpolation**, optional: if True, interpolate values.

        :Return: A :class:`~configobj.ConfigObj` instance
        """
        if interpolation is None or interpolation is True: interpolation = self._interpolation
        cfg = ConfigObj(interpolation=interpolation, configspec=self._configspec,
            encoding=self._encoding)
        cfg.validate(self._validator, copy=True)
        if nocomments:
            cfg.walk(_walker_remove_all_comments_, call_on_sections=int(nocomments)==2)
        elif self._configspec:
            cfg.inline_comments = self._configspec.inline_comments
        return cfg

    def reset(self, cfgfile='config.cfg', backup=True, nocomments=True, verbose=True):
        """Reset a config file  to default values

        :Params:

            - **cfgfile**, optional: The configuration file to reset.
            - **backup**, optional: Backup the old config file.
            - **nocomments**, optional: Do not include comment in config file.

        :Return: A :class:`~configobj.ConfigObj` instance

        :See also: :meth:`defaults`
        """

        # Load defaults
        cfg = self.defaults(nocomments=nocomments, interpolation=False)

        # Remove old file
        if os.path.exists(cfgfile):
            if backup:
                shutil.copy(cfgfile, cfgfile+'.bak')
            os.remove(cfgfile)
        else:
            backup = False

        # Write to new one
        cfg.filename = cfgfile
        cfg.write()
        if verbose:
            print 'Created default config file: %s'%cfgfile
            if backup:
                print "Backuped old config file to: %s.bak"%cfgfile
        return cfg

    def load(self, cfgfile='config.cfg',
        validate='fix', geterr=False, patch=None, force=True, **kwpatch):
        """Get a :class:`~configobj.ConfigObj` instance loaded from a file

        :Params:

            - **cfgfile**, optional: config file

                - a config file name
                - a :class:`~configobj.ConfigObj` instance
                - ``None``: defaults to ``"config.cfg"``

            - **validate**, optional: Type of validation

                - ``False``: no validation
                - ``"fix"``: validation fixes and reports errors
                - ``"report"``: validation reports errors
                - ``"raise"``: validation raises errors

            - **geterr**, optional: Return validation results as well
            - **force**, optional: Force re-instantiation of ``cfgfile`` when it is already
              a :class:`ConfigObj` instance.

        :Return: Depends on ``geterr``

            - if ``True``: ``(cfg, err)`` where is the result of :meth:`~configobj.ConfigObj.validate`
            - else: ``cfg`` (:class:`~configobj.ConfigObj` instance)
        """

        # Load the config
        if cfgfile is not None and isinstance(cfgfile, basestring) and not os.path.exists(cfgfile):
            cfgfile = None

        # Instantiate / Copy
        if not isinstance(cfgfile, ConfigObj) or force:
            cfg = ConfigObj(cfgfile, interpolation=self._interpolation,
                configspec=self._configspec, encoding=self._encoding)
        else:
            cfg = cfgfile


        # Patch
        if kwpatch and not patch:
            patch = {}
        if patch is not None:
            if kwpatch: # Patch the patch!
                patch = self.patch(patch, kwpatch, validate=False)
            self.patch(cfg, patch, validate=False)

        # Validation
        if validate and self._configspec:

            # Validation itself
            err = cfg.validate(self._validator, preserve_errors=True)

            # Defaults for fixing
            if validate in ['fix', 'report']:
                defaults = self.defaults()

            # Loop on errors
            if isinstance(err, dict):

                for sections, key, error in flatten_errors(cfg, err):

                    # Format explicit message
                    if len(sections):
                        section = '['+']['.join(sections)+'] '
                    else:
                        section = ''
                    msg = 'Config value error: %s%s: %s'%(section, key, getattr(error, 'message', error))

                    # Check what to do
                    if validate in ['fix', 'report']:


                        # Report
                        sec = cfg
                        secd = defaults
                        for subsec in sections:
                            sec = sec[subsec]
                            if subsec not in secd:
                                secd = None
                                break
                            secd = secd[subsec]
                        if secd is None:
#                            sys.stderr.write(msg+"\nCan't set it to default value because"
#                                " none available")
                            continue

                        # Loop on keys
                        keys = secd.keys() if key is None else [key]
                        for k in keys:

                            vdef = secd.get(k, None)
                            sys.stderr.write(msg+'\nSetting it to default value: %s\n'%vdef)#,ValidationWarning)

                            # Reset to default
                            if validate=='fix':
                                sec[k] = vdef

                    else: # Raise an error
                        raise ValidateError, msg

            if geterr:
                return cfg, err

        return cfg


    def patch(self, cfg, cfgpatch, validate=False):
        """Replace config values of ``cfg`` by those of ``cfgpatch``

        :Params:

            - **cfg**: A :class:`~configobj.ConfigObj` instance, a config file or
              a dictionary, that must be patched.
            - **cfgpatch**: A :class:`~configobj.ConfigObj` instance, a config file or
              a dictionary, used for patching.
            - **validate**, optional: If ``True``, validate configs if they have a valid config spec.
        """
        if not isinstance(cfg, ConfigObj):
            cfg = ConfigObj(cfg, configspec=self._configspec, encoding=self._encoding)#, interpolation=False)
        #else:
            #cfg.interpolation = False
        if not isinstance(cfgpatch, ConfigObj):
            cfgpatch = ConfigObj(cfgpatch, configspec=self._configspec, interpolation=False,
                encoding=self._encoding)
        else:
            cfgpatch.interpolation = False
        cfgpatch.walk(_walker_patch_, cfg=cfg)
        if validate and cfg.configspec is not None:
            cfg.validate(self.validator())
        #cfg.interpolation = 'template'
        return cfg

    def arg_parse(
            self, parser=None, exc=[], parse=True, args=None,
            getparser=False, getargs=False, cfgfile='config.cfg',
            patch=None, cfgfileopt='--cfgfile', cfgfilepatch='before',
            nested=None):
        """Options (:mod:`argparse`) and config mixer.

            1. Creates command-line options from config defaults
            2. Parse command-line argument and create a configuration patch

        For instance, the following config define the commandline option
        ``--section1-my-section2-my-key`` with ``value`` as a default value, s
        tored in a special group of options with a short name and a long description::

            [section1] # Short name : long description of the group
                [[my_section2]]
                    my_key=value

        .. warning:: Section and option names must not contain any space-like character !

        .. note::

            If you want to prevent conflict of options, don't use ``"_"`` in
            section and option names.


        :Params:

            - **parser**: optional, a default one is created if not given. This can be:
                - a :class:`OptionParser` instance
                - a :class:`dict` with keyword arguments for the one to be created
            - **exc**, optional: List of keys to be excluded from parsing.
            - **parse**, optional: If ``True``, parse commande line options and arguments
            - **args**, optional: List of arguments to parse instead of default sys.argv[1:]
            - **getparser**: allow getting the parser in addition to the config if parse=True
            - **getargs**: allow getting the parsed arguments in addition to the config if parse=True
            - **patch**, optional: used if parse is True.
                Can take the following values:
                - a :class:`bool` value indicating wheter to apply defaults on the returned config
                  before applying the command line config
                - a :class:`ConfigObj` instance to apply on the returned config
                  before applying the command line config
            - **cfgfileopt**: optional, if present a config file option will be added.
                Can be a :class:`string` or couple of strings to use as the option short and/or long name
            - **cfgfilepatch**: specify if the returned config must be patched if a config file command
                line option is provided and when to patch it. Can take the following values:
                - True or 'before': the config file would be used before command line options
                - 'after': the config file would be used after command line options
                - Any False like value: the config file would not be used
            - **nested**: Name of a section whose defines the configuration.
              It must be used when the configuration in nested in more general configuration.

        :Return:
            - the :class:`OptionParser` if parse is False
            - the resulting :class:`~configobj.ConfigObj` if parse is True and getparser is not True
            - the resulting :class:`~configobj.ConfigObj` and the :class:`OptionParser` if both parse and getparser are True
        """

        # Prepare the option parser
        if parser is None:
            parser = ArgumentParser(add_help=False)
        elif isinstance(parser, dict):
            parser['add_help'] = False
            parser = ArgumentParser(**parser)

        # Add short and long helps
        old_conflict_handler = parser._optionals.conflict_handler
        parser._optionals.conflict_handler = 'resolve'
        parser.add_argument('-h','--help', action=AP_ShortHelpAction,
            help='show a reduced help and exit')
        parser.add_argument('--long-help', action=_HelpAction,
            help='show an extended help and exit')
        parser._optionals.conflict_handler = old_conflict_handler

        # Add the cfgfile option (configurable)
        if cfgfileopt:
            if isinstance(cfgfileopt, basestring):
                if not cfgfileopt.startswith('-'):
                    if len(cfgfileopt)==1:
                        cfgfileopt = '-'+cfgfileopt
                    else:
                        cfgfileopt = '--'+cfgfileopt
                cfgfileopt = (cfgfileopt,)
            parser.add_argument(*cfgfileopt,
                dest="cfgfile", help='user configuration file that overrides defauts'
                    '[default: "%(default)s"]', default=cfgfile)

        # Default config
        defaults = self.defaults()

        # Create global group of options from defaults
        # - inits
        re_match_initcom = re.compile(r'#\s*-\*-\s*coding\s*:\s*\S+\s*-\*-\s*').match
        if len(defaults.initial_comment)==0 or re_match_initcom(defaults.initial_comment[0]) is None:
            desc = ['global configuration options']
        else:
            re_match_initcom(defaults.initial_comment[0]), defaults.initial_comment[0]
            icom = int(re_match_initcom(defaults.initial_comment[0]) is not None)
            if len(defaults.initial_comment)>icom:
                desc = defaults.initial_comment[icom].strip('# ').split(':', 1)
        group = parser.add_argument_group(*desc)
        # - global options
        for key in defaults.scalars:
            if key not in exc:
                _walker_argcfg_setarg_(defaults, key, group=group, exc=exc, nested=nested,
                    boolean_false=self._boolean_false)
#                group.add_argument('--'+_cfg2optname_(key, nested), help=_shelp_(defaults, key))
            else:
                pass

        # Create secondary option groups from defaults
        for key in defaults.sections:
            desc = [key]
            comment = defaults.inline_comments[key] # FIXME: always empty!
            if comment is not None:
                desc = comment.strip('# ')
                if ':' in desc:
                    desc = desc.split(':', 1)
                else:
                    desc = [key.lower(), desc]
            section = defaults[key]
            group = parser.add_argument_group(*desc)
            defaults[key].walk(_walker_argcfg_setarg_, raise_errors=True,
                call_on_sections=False, group=group, exc=exc, nested=nested,
                boolean_false=self._boolean_false)

        # Now create a configuration instance from passed options
        if parse:

            # Which args ?
            if args is None: args = sys.argv[1:]

            # Parse
            options = parser.parse_args(list(args))

            # Create a configuration to feed
            cfg = ConfigObj(interpolation=self._interpolation, encoding=self._encoding)

            # Initial config from defaults or the one supplied
            if patch:
                self.patch(cfg, patch if isinstance(patch, ConfigObj) else defaults)
            if cfgfilepatch:
                if (isinstance(cfgfilepatch, basestring)
                and cfgfilepatch.strip().lower().startswith('a')):
                    cfgfilepatch = 'after'
                else:
                    cfgfilepatch = 'before'

            # Feed config with cfgfile before command line options
            if cfgfilepatch == 'before' and getattr(options, 'cfgfile', None):
                self.patch(cfg, self.load(options.cfgfile))

            # Feed config with command line options
            defaults.walk(_walker_argcfg_setcfg_, raise_errors=True,
                call_on_sections=False, cfg=cfg, options=options, nested=nested)

            # Feed config with cfgfile after command line options
            if cfgfilepatch == 'after' and getattr(options, 'cfgfile', None):
                self.patch(cfg, self.load(options.cfgfile))

            if not getparser and not getargs: return cfg
            out = cfg,
            if getparser: out += parser,
            if getargs: out += options,
            return out

        else:
            return parser

    def arg_patch(self, parser, exc=[], cfgfileopt='cfgfile'):
        """Call to :meth:`arg_parse` with an automatic patching of current configuration

        :Return: ``cfg, args``
        """

        # Create a patch configuration from commandline arguments
        cfgpatch, args = self.arg_parse(parser, exc=exc, cfgfileopt=cfgfileopt, getargs=True)

        #  Load personal config file and default values
        cfg = self.load(args.cfgfile)

        #  Patch it with commandeline options
        self.patch(cfg, cfgpatch)

        return cfg, args


    def opt_parse(
            self, parser=None, exc=[], parse=True, args=None,
            getparser=None, cfgfile='config.cfg',
            patch=None, cfgfileopt='--cfgfile', cfgfilepatch='before',
            nested=None):
        """Options (:mod:`optparse`) and config mixer.

            1. Creates command-line options from config defaults
            2. Parse command-line argument and create a configuration patch

        For instance, the following config define the commandline option
        ``--section1-my-section2-my-key`` with ``value`` as a default value, s
        tored in a special group of options with a short name and a long description::

            [section1] # Short name : long description of the group
                [[my_section2]]
                    my_key=value

        .. warning:: Section and option names must not contain any space-like character !

        .. note::

            If you want to prevent conflict of options, don't use ``"_"`` in
            section and option names.


        :Params:

            - **parser**: optional, a default one is created if not given. This can be:
                - a :class:`OptionParser` instance
                - a :class:`dict` with keyword arguments for the one to be created
            - **exc**, optional: List of keys to be excluded from parsing.
            - **parse**, optional: If ``True``, parse commande line options and arguments
            - **args**, optional: List of arguments to parse instead of default sys.argv[1:]
            - **getparser**: allow getting the parser in addition to the config if parse=True
            - **patch**, optional: used if parse is True.
                Can take the following values:
                - a :class:`bool` value indicating wheter to apply defaults on the returned config
                  before applying the command line config
                - a :class:`ConfigObj` instance to apply on the returned config
                  before applying the command line config
            - **cfgfileopt**: optional, if present a config file option will be added.
                Can be a :class:`string` or couple of strings to use as the option short and/or long name
            - **cfgfilepatch**: specify if the returned config must be patched if a config file command
                line option is provided and when to patch it. Can take the following values:
                - True or 'before': the config file would be used before command line options
                - 'after': the config file would be used after command line options
                - Any False like value: the config file would not be used

        :Return:
            - the :class:`OptionParser` if parse is False
            - the resulting :class:`~configobj.ConfigObj` if parse is True and getparser is not True
            - the resulting :class:`~configobj.ConfigObj` and the :class:`OptionParser` if both parse and getparser are True
        """

        # Prepare the option parser
        if parser is None:
            parser = OptionParser()
        elif isinstance(parser, dict):
            parser = OptionParser(**parser)

        # Add (or override!) option types checkers
        # Define a new Option class
        class option_class(Option):
            TYPES = tuple(list(parser.option_class.TYPES) + list(self._validator.functions.keys()))
            TYPE_CHECKER = parser.option_class.TYPE_CHECKER.copy()
        # Define the wrapping function for optparse option types which also handle list values
        def wrap_option_type_checker(func, islist):
            def wrapper_option_type_checker(opt, name, value):
                if islist: # Use configobj list parser
                    value,comment = ConfigObj(list_values=True, interpolation=False,
                        encoding=self._encoding)._handle_value(value)
                    return func(value)
                else:
                    return func(value)
            wrapper_option_type_checker.__name__ += '-'+func.__name__
            return wrapper_option_type_checker
        # Setup type checkers  into the new Option class
        for name,func in self._validator.functions.items():
            if name in option_class.TYPE_CHECKER:
                pass # warn('Overriding Option type checker %s'%(name))
            islist = _validator_specs_.get(name, {}).get('iterable', None)
            option_class.TYPE_CHECKER[name] = wrap_option_type_checker(func, islist)
        # Replace the parser Option class
        parser.option_class = option_class

        # Add the cfgfile option (configurable)
        if cfgfileopt:
            if isinstance(cfgfileopt, basestring):
                cfgfileopt = (cfgfileopt,)
            parser.add_option(*cfgfileopt, action='store', type="string",
                dest="cfgfile", help='Configuration file [default: "%s"]'%cfgfile, default=cfgfile)

        # Default config
        defaults = self.defaults()

        # Create global group of options from defaults
        # - inits
        re_match_initcom = re.compile(r'#\s*-\*-\s*coding\s*:\s*\S+\s*-\*-\s*').match
        if len(defaults.initial_comment)==0 or re_match_initcom(defaults.initial_comment[0]) is None:
            desc = ['Global configuration options']
        else:
            re_match_initcom(defaults.initial_comment[0]), defaults.initial_comment[0]
            icom = int(re_match_initcom(defaults.initial_comment[0]) is not None)
            if len(defaults.initial_comment)>icom:
                desc = defaults.initial_comment[icom].strip('# ').split(':', 1)
        group = OptionGroup(parser, *desc)
        # - global options
        for key in defaults.scalars:
            if key not in exc:
                _walker_optcfg_setopt_(defaults, key, group=group, exc=exc, nested=nested,
                    boolean_false=self._boolean_false)
#                group.add_option('--'+_cfg2optname_(key, nested), action='store', type="string",
#                    dest=key, help=_shelp_(defaults, key))
            else:
                pass
        # - add to parser
        parser.add_option_group(group)

        # Create secondary option groups from defaults
        for key in defaults.sections:
            desc = ['Undocumented section']
            comment = defaults.inline_comments[key]
            if comment is not None:
                desc = comment.strip('# ').split(':', 1)
            section = defaults[key]
            group = OptionGroup(parser, *desc)
            defaults[key].walk(_walker_optcfg_setopt_, raise_errors=True,
                call_on_sections=False, group=group, exc=exc, nested=nested,
                boolean_false=self._boolean_false)
            parser.add_option_group(group)

        # Now create a configuration instance from passed options
        if parse:

            if args is None: args = sys.argv[1:]
            (options, args) = parser.parse_args(list(args))

            # Intercept helps
            if getattr(options, 'long_help', None):
                parser.print_help()
                sys.exit()
            elif getattr(options, 'help', None):
                print_short_help(parser)
                sys.exit()

            # Create a configuration to feed
            cfg = ConfigObj(interpolation=self._interpolation, encoding=self._encoding)

            # Initial config from defaults or the one supplied
            if patch:
                self.patch(cfg, patch if isinstance(patch, ConfigObj) else defaults)
            if cfgfilepatch:
                if (isinstance(cfgfilepatch, basestring)
                and cfgfilepatch.strip().lower().startswith('a')):
                    cfgfilepatch = 'after'
                else:
                    cfgfilepatch = 'before'

            # Feed config with cfgfile before command line options
            if cfgfilepatch == 'before' and getattr(options, 'cfgfile', None):
                self.patch(cfg, self.load(options.cfgfile))

            # Feed config with command line options
            defaults.walk(_walker_optcfg_setcfg_, raise_errors=True,
                call_on_sections=False, cfg=cfg, options=options, nested=nested)

            # Feed config with cfgfile after command line options
            if cfgfilepatch == 'after' and getattr(options, 'cfgfile', None):
                self.patch(cfg, self.load(options.cfgfile))
            return (cfg, parser) if getparser else cfg
        else:
            return parser

    def opt_patch(self, parser, exc=[], cfgfileopt='cfgfile'):
        """Call to :meth:`arg_parse` with an automatic patching of current configuration

        :Return: ``cfg, args``
        """

        # Create a patch configuration from commandline arguments
        cfgpatch = self.opt_parse(parser, exc=exc, cfgfileopt=cfgfileopt)

        #  Load personal config file and default values
        cfg = self.load(parser.values.cfgfile)

        #  Patch it with commandeline options
        self.patch(cfg, cfgpatch)

        return cfg





    def opt_long_help(self, rst=True, usage="Usage: [prog] [options] ...",
        description="Long help based on config specs"):
        """Get the generic long help from config specs

        :Params:

            - **rst**, optional: Reformat output in rst.
        """

        # Standard options
        parser = OptionParser(usage=usage, description=description, add_help_option=False )
        parser.add_option('-h','--help', action='store_true', dest="help",
            help='show a reduced help')
        parser.add_option('--long-help', action='store_true', dest="long_help",
            help='show an extended help')

        # Configuration options
        self.opt_parse(parser, parse=False)
        if rst:
            formatter = IndentedHelpFormatter(width=1000)#max_help_position=0)
            shelp = parser.format_option_help(formatter).encode(sys.getdefaultencoding(), 'replace')
        else:
            shelp = parser.format_help().encode(sys.getdefaultencoding(), 'replace')

        # Convert to rst
        if rst: shelp = opt2rst(shelp)

        return shelp

    def arg_long_help(self, rst=True, usage=None,
        description="Long help based on config specs"):
        """Get the generic long help from config specs

        :Params:

            - **rst**, optional: Reformat output in rst.
        """

        # Standard options
        parser = ArgumentParser(usage=usage, description=description, add_help=False )
        parser.add_argument('-h','--help', action='store_true', help='show a reduced help')
        parser.add_argument('--long-help', action='store_true', help='show an extended help')

        # Configuration options
        self.arg_parse(parser, parse=False)
        if rst:
            formatter = ArgHelpFormatter(max_help_position=0)
            for action_group in parser._action_groups:
                formatter.start_section(action_group.title)
                formatter.add_text(action_group.description)
                formatter.add_arguments(action_group._group_actions)
                formatter.end_section()
            shelp = formatter.format_help()
        else:
            shelp = parser.format_help()

        # Encoding
        shelp = shelp.encode(sys.getdefaultencoding(), 'replace')

        # Convert to rst
        if rst: shelp = opt2rst(shelp)

        return shelp

    def get_rst(self):
        """Convert the default configuration to rst declarations with :func:`cfg2rst`"""
        return cfg2rst(self)


def cfgargparse(cfgspecfile, parser, cfgfileopt='cfgfile', cfgfile='config.cfg',
    exc=[], **kwargs):
    """Merge configuration and commandline arguments

    :Params:

        - **cfgspecfile**: Config specification file (.ini).
        - **parser**: :class:`~argpase.ArgumentParser` instance.
        - **cfgfileopt**, optional: Name of the option used to specify the
          user config file. Example: ``'cfgfile'`` creates the option
          ``--cfgfile=<config file>``.
        - **cfgfile**, optional: Default name for the loaded config file.
        - **exc**, optional: Config option name that must not be used to generated
          a commandline option.
        - Extra params are passed to :class:`ConfigManager` initialization.

    :Return: A :class:`ConfigObj` object

    :Tasks:

        1. Initialize a default configuration (:class:`ConfigManager`)
           from the specification file given by ``cfgspecfile``.
        2. Generate associated commandline options.
        3. Load a user configuration file (specified with the option
           whose name is given by ``cfgfileopt``).
        4. Patch this configuration with user supplied options retrieved
           using the :class:`~argpase.ArgumentParser` parser ``parser``.

        Technically it combines :class:`ConfigManager` and :meth:`ConfigManager.arg_parse`
    """
    return ConfigManager(cfgspecfile, **kwargs).arg_parse(
        parser, cfgfileopt=cfgfileopt, exc=exc, cfgfile=cfgfile, getargs=True)

def cfgoptparse(cfgspecfile, parser, cfgfileopt='cfgfile', cfgfile='config.cfg',
    exc=[], **kwargs):
    """Merge configuration and commandline arguments

    :Params:

        - **cfgspecfile**: Config specification file (.ini).
        - **parser**: :class:`~argpase.ArgumentParser` instance.
        - **cfgfileopt**, optional: Name of the option used to specify the
          user config file. Example: ``'cfgfile'`` creates the option
          ``--cfgfile=<config file>``.
        - **cfgfile**, optional: Default name for the loaded config file.
        - **exc**, optional: Config option name that must not be used to generated
          a commandline option.
        - Extra params are passed to :class:`ConfigManager` initialization.

    :Tasks:

        1. Initialize a default configuration (:class:`ConfigManager`)
           from the specification file given by ``cfgspecfile``.
        2. Generate associated commandline options.
        3. Load a user configuration file (specified with the option
           whose name is given by ``cfgfileopt``).
        4. Patch this configuration with user supplied options retrieved
           using the :class:`~optpase.OptionParser` parser ``parser``.

        Technically it combines :class:`ConfigManager` and :meth:`ConfigManager.opt_patch`
    """
    return ConfigManager(cfgspecfile, **kwargs).opt_parse(
        parser, cfgfileopt=cfgfileopt, exc=exc, cfgfile=cfgfile, getparser=True)


def opt2rst(shelp, prog=None, secfmt=':%(secname)s:', descname='Description'):
    """Convert options displayed in an help string to rst declarations of options

    This is useful for autodocumenting executable python scripts that show a formatted help.

    :Params:

        - **shelp**: Help string showing options (results from :option:``--help``).
        - **prog**, optional: Program name, otherwise guess it from usage.

    :Output: String converted to rst format (with :rst:dir:`cmdoption` directives).

    """
    rhelp = []
    multiline = False
    s_param = r'(?:\{[\w,]+\}|\w+)'
    s_sopt = r'(?:-\w+(?: %(s_param)s)?)'%locals() # short option (-t)
    s_lopt = r'(?:--[\w\-]+(?:[= ]+%(s_param)s)?)'%locals() # long option (--toto)
    s_optsep = r'(?:, +)' # option separator
    s_desc = r'(?:  (.+))'
    s_tot = r'^  (?:  )?((?:%(s_sopt)s|%(s_lopt)s)(?:%(s_optsep)s(?:%(s_sopt)s|%(s_lopt)s))*)%(s_desc)s?$'%locals()
#    s_tot = r'^  (%(s_sopt)s|%(s_lopt)s)|%(s_sopt)s%(s_optsep)s%(s_lopt)s)%(s_desc)s?$'%locals()
    re_opt = re.compile(s_tot).match
    re_sec = re.compile(r'^(?:  )?([\w\s]+):(?: (.+))?$').match
    secname = None
    for line in shelp.splitlines():

        # Sections
        m = re_sec(line)
        if m and not line.lower().endswith('ex:'):

            secname = m.group(1).title().strip()

            # Usage
            if secname=='Usage' and m.group(2) is not None:
                usage = m.group(2).strip()
                if prog is None:
                    prog = os.path.basename(usage.split()[0])
                rhelp.append('.. program:: %s\n'%prog)
                rhelp.extend([secfmt%locals(), "\n\t.. code-block:: bash\n\n\t\t%s"%usage])
                multiline = True
            else:
                rhelp.extend([secfmt%locals(), ''])
                if m.group(2) is not None:
                    rhelp.extend(['', '\t'+m.group(2)])
                    multiline = True
            continue

        # Options and other lines
        m = re_opt(line)
        if m:

            rhelp.extend(['','\t.. cmdoption:: '+m.group(1), ''])
            multiline = True
            if m.group(2) is not None:
                rhelp.append('\t\t'+m.group(2).strip())

        elif secname and secname.lower()=='positional arguments' and line.startswith(' '*2):

            sline = line.split()
            rhelp.extend(['','\t.. cmdoption:: '+sline[0], ''])
            multiline = True
            if len(sline)>1 is not None:
                rhelp.append('\t\t'+' '.join(sline[1:]))

        elif multiline and len(line.strip()) and line.startswith(' '*3):

            indent = '\t\t'
            if secname=='Usage':
                indent += '\t'
            rhelp.append(indent+line.strip())
        #elif secname==descname:
        #    rhelp.append('\t'+line)

        else:

            indent = ''
            if secname==descname:
                indent += '\t'
            rhelp.append(indent+line)
            multiline = False
            if secname=='Usage':
                secname = descname
                rhelp.extend([secfmt%locals(), ''])

    return '\n'.join(rhelp)

def _opt2cfgname_(name, nested):
    cfgkey = name.replace('-', '_')
    if nested and cfgkey.startswith(nested+'_'):
        cfgkey = cfgkey[len(nested+'_'):]
    return cfgkey

_re_cfg2optname_sub = re.compile('[_\s]').sub
def _cfg2optname_(name, nested=None):
    optkey = _re_cfg2optname_sub('-', name)
    if nested: optkey = nested+'-'+optkey
    return optkey.lower()

class _attdict_(dict):
    def __getattr__(self, name):
        if name in self.__dict__: return object.__getattribute__(self, name)
        else: return self[name]

def getspec(spec, validator=None):
        '''
        Get an option specification.

        :Params:
            - **spec**: the specification string
            - **validator**: (optional) the validator to use

        :Return: A dict with keys:
            - **funcname**: the validation function name
            - **args**:  the positionnal arguments
            - **kwargs**: the named arguments
            - **default**: the default value
            - **iterable**: if the value is list-like
            - **func**: the validation function
            - **opttype**: the function used with :mod:`optparse`
            - **argtype**: the function used with :mod:`argparse`

        Read access to these keys can also be done as attribute of the returned dict
        (d.funcname == d['funcname'], ...)

        For example, a specification file containing:

        [section]
            option = integer(default=0, min=-10, max=10)

        Would return: (integer, [], {'min': '-10', 'max': '10'}, 0)

        This can be usefull when you added extraneous named arguments into your
        specification file for your own use.

        '''
        if not validator:
            validator = Validator()
            validator.functions.update(_validator_functions_)
        funcname, args, kwargs, default = validator._parse_with_caching(spec)
        spec = _validator_specs_.get(funcname, dict(func=None, iterable=None, opttype=None, argtype=None)).copy()
        spec.update(dict(funcname=funcname, args=args, kwargs=kwargs, default=default))
        return _attdict_(spec)

def _walker_remove_all_comments_(sec, key):
    sec.comments[key] = ''
    sec.inline_comments[key] = ''

def _walker_remove_comments_(sec, key):
    sec.comments[key] = ''

def _walker_remove_inline_comments_(sec, key):
    sec.inline_comments[key] = ''

def _walker_unchanged_options_(sec, key):
    if not sec.configspec: return
    spec = getspec(sec.configspec.get(key, ''))
    return spec.default == sec[key]

def remove_defaults(cfg):
    defaults = cfg.walk(_walker_unchanged_options_, call_on_sections=False)
    def remove(c, d):
        for k,v in d.items():
            if isinstance(v, dict):
                remove(c[k], v)
            elif v:
                c.pop(k)
    remove(cfg, defaults)

def _walker_optcfg_setcfg_(sec, key, cfg=None, options=None, nested=None):
    """Walker to set config values"""
    # Find genealogy
    parents = _parent_list_(sec, names=False)
    cfgkey = '_'.join([p.name.strip('_') for p in parents]+[key])
    for option, value in options.__dict__.items():

        # Option not set
        if value is None: continue

        # Option matches key?
        if _opt2cfgname_(option, nested) != cfgkey.lower(): continue

        # Check or create cfg genealogy
        s = cfg
        for p in parents:
            if not s.has_key(p.name):
                s[p.name] = {}
            s = s[p.name]
        s[key] = value

def _walker_optcfg_setopt_(sec, key, group=None, exc=None, nested=None, boolean_false=True):
    """Walker to set options"""
    # Find option key and output var name
    key = key.strip('_')
    pp = _parent_list_(sec)
    varname = '_'.join(pp+[key])
    optkey = _cfg2optname_(varname, nested)

    # Check exceptions
    if key in exc: return

    # Add option to group
    spec = _validator_specs_.get(sec.configspec[key].split('(', 1)[0], {})
    type = spec.get('opttype', 'string')
    if boolean_false and type=='boolean':
        default = sec[key]
        action = 'store_false' if default else 'store_true'
        type = None
    else:
        action = 'store'
        default = None
    group.add_option(
        '--'+optkey,
            #action='append' if spec.get('iterable', None) else 'store',
            action=action,
            type=type,
            dest=varname,
            help=_shelp_(sec, key),
            default=default,
    )

def _walker_argcfg_setcfg_(sec, key, cfg=None, options=None, nested=None):
    """Walker to set config values"""
    # Find genealogy
    parents = _parent_list_(sec, names=False)
    cfgkey = '_'.join([p.name.strip('_') for p in parents]+[key])
    for option, value in options._get_kwargs():

        # Option not set
        if value is None: continue

        # Option matches key?
        if _opt2cfgname_(option, nested) != cfgkey.lower(): continue

        # Check or create cfg genealogy
        s = cfg
        for p in parents:
            if not s.has_key(p.name):
                s[p.name] = {}
            s = s[p.name]
        s[key] = value

def _walker_argcfg_setarg_(sec, key, group=None, exc=None, nested=None, encoding=None,
    boolean_false=True):
    """Walker to set options"""
    # Find option key and output var name
    key = key.strip('_')
    pp = _parent_list_(sec)
    varname = '_'.join(pp+[key])
    optkey = _cfg2optname_(varname, nested)

    # Check exceptions
    if key in exc: return
    spec = _validator_specs_.get(sec.configspec[key].split('(', 1)[0], {})

    # Define the wrapping function for argparse argument types which also handle list values
    def wrap_argparse_type(func, islist):
        def wrapper_argparse_type(value):
            if islist: # Use configobj list parser
                value,comment = ConfigObj(list_values=True, interpolation=False,
                    encoding=encoding)._handle_value(value)
                return func(value)
            else:
                return func(value)
        wrapper_argparse_type.__name__ += '-'+func.__name__
        return wrapper_argparse_type

    # Add argument to group
    type = wrap_argparse_type(spec.get('func', lambda s:s), spec.get('iterable', None))
    kw = {}
    if boolean_false and spec.get('opttype', '')=='boolean':
        default = sec[key]
        action = 'store_true' if default is False else 'store_false'
    else:
        action = 'store'
        kw['type'] = type
    group.add_argument(
        '--'+optkey,
            action=action,
            help=_shelp_(sec, key),
            **kw
    )

def _shelp_(sec, key, format='%(shelp)s [default: %(default)r]', mode='auto',
    undoc='Undocumented', adddot=True):
    """Get help string

    :Params:

        - **mode**:

            - inline: inline comment only,
            - above: above comments only,
            - merge: merge inline and above comments,
            - auto: if one is empty use the other one, else use inline

    """
    # filter
    def strip(c):
        return c.strip().strip('#').strip()
    abcoms = map(strip, filter(lambda c: c is not None, sec.comments[key]))
    incoms = map(strip, filter(lambda c: c is not None, [sec.inline_comments[key]]))


    # Merge comments above item and its inline comment
    if mode=='merge':
        comments = abcoms+incoms
    elif mode=='auto':
        if not incoms:
            comments = abcoms
        else:
            comments = incoms
    elif mode=='above':
        comments = abcoms
    else:
        comments = incoms

    # Force comments to end with a dot '.'
    if adddot:
        comments = [c.endswith('.') and c or '%s.'%c for c in comments]
    shelp = '\n'.join(comments)

    # If no comments
    if not shelp: shelp = undoc
    default = _sdefault_(sec, key)
    return format%locals()

def _sdefault_(sec, key):
    """Get default string"""
    default = sec.get(key, None)
    if isinstance(default, (list,tuple)):
        default = ','.join(map(str, default))
    else:
        default = str(default).strip('()')
    return default.replace('%', '%%')

def _parent_list_(sec, names=True):
    parents = []
    while sec.name is not None:
        parents.insert(0, names and sec.name.strip('_') or sec)
        sec = sec.parent
    return parents

def _walker_patch_(patch_sec, patch_key, cfg):
    """Walk through the patch to apply it"""
    sec = cfg
    for key in _parent_list_(patch_sec):
        if not sec.has_key(key):
            sec[key] = {}
        sec = sec[key]
    try:
        patch_sec[patch_key] = type(sec[patch_key])(patch_sec[patch_key])
    except:
        pass
    sec[patch_key] = patch_sec[patch_key]

def _walker_cfg2rst_(cfg, key, lines):

    # Name
    secnames = get_secnames(cfg)
    if key in cfg.sections:
        secnames.append(key)
        name = '[%s]'%(']['.join(secnames), )
        conftype = 'confsec'
    else:
        if secnames:
            name = '[%s] %s'%(']['.join(secnames), key)
        else:
            name = key
        conftype = 'confopt'

    # Description
    desc = cfg.inline_comments.get(key)
    if desc is None: desc = ''
    desc = desc.strip('#').strip()

    # Formatting
    desc = redent(desc, 1)
    text = '.. %(conftype)s:: %(name)s\n\n%(desc)s\n'%locals()
    text = redent(text, cfg.depth)
    lines.append(text)


def _walker_set_boolean_false_by_default_(sec, key, validator=None):
    if validator is None: return
    check = sec[key]
    fun_name, fun_args, fun_kwargs, default = validator._parse_with_caching(check)
    if fun_name=='boolean' and default is None:
        if fun_args or fun_kwargs:
            check = check[:-1]+', default=False)'
        else:
            check = check+'(default=False)'
        sec[key] = check

def get_secnames(cfg):
    """Get section names as list from top to bottom ['sec0','sec1',...]"""
    if cfg.depth==0: return []
    secnames = [cfg.name]
    for i in xrange(cfg.depth-1):
        cfg = cfg.parent
        secnames.append(cfg.name)
    return secnames[::-1]

def redent(text, n=1, indent='    '):
    lines = text.split('\n')
    lines = [(n*indent+line) for line in lines]
    return '\n'.join(lines)


def cfg2rst(cfg):
    """Convert a configuration to rst format

    Configuration sections are declared with the rst directive
    "confsec" and options are declared with the rst directive
    "confopt".

    For instance:

    .. code-block:: ini

        a=1 # desc a
        [s1] # desc s1
            b=2  # desc b
            [[s2]] # desc s2
                c=$a-$b # desd c
                [sec1] # section 1

    is converted to:

    .. code-block:: rst

        .. confopt:: a

            desc a

        .. confsec:: [s1]

            desc s1

            .. confopt:: [s1] b

                desc b

            .. confsec:: [s1][s2]

                desc s2

                .. confopt:: [s1][s2] c

                    desd c

    Then one can reference an option with for example ``:confopt:`[s1][s2]c`.

    :Params:

        - **cfg**: :class:`ConfigObj` or :class:`ConfigManager` instance.
          In the case of a :class:`ConfigManager`, the :meth:`~ConfigManager.defaults`
          are used.

    :Return: rst string
    """
    out = ''
    if isinstance(cfg, ConfigManager):
        cfg = cfg.defaults()
    lines = []
    cfg.walk(_walker_cfg2rst_, call_on_sections=True, lines=lines)
    return '\n'.join(lines)


def print_short_help(parser, formatter=None):
    """Print all help of an :class:`~optparse.OptionParser` instance but those of groups."""
    if isinstance(parser, OptionParser):
        # - top
        if formatter is None:
            formatter = parser.formatter
        result = []
        if parser.usage:
            result.append(parser.get_usage() + "\n")
        if parser.description:
            result.append(parser.format_description(formatter) + "\n")
        # - basic options skipping groups
        formatter.store_option_strings(parser)
        result.append(formatter.format_heading(_("Options")))
        formatter.indent()
        if parser.option_list:
            result.append(OptionContainer.format_option_help(parser, formatter))
            result.append("\n")
        formatter.dedent()
        del result[-1]
        result.append(parser.format_epilog(formatter))
        print "".join(result)

    elif isinstance(parser, ArgumentParser):

        if formatter is None:
            formatter = parser._get_formatter()

        # usage
        formatter.add_usage(parser.usage, parser._actions,
                            parser._mutually_exclusive_groups)

        # description
        formatter.add_text(parser.description)

        # positionals, optionals and user-defined groups
        for action_group in parser._action_groups:
            if action_group.title in ['positional arguments', 'optional arguments']:
                formatter.start_section(action_group.title)
                formatter.add_text(action_group.description)
                formatter.add_arguments(action_group._group_actions)
                formatter.end_section()

        # epilog
        formatter.add_text(parser.epilog)

        # determine help from format above
        parser._print_message(formatter.format_help(), sys.stdout)


class AP_ShortHelpAction(_HelpAction):

    def __call__(self, parser, namespace, values, option_string=None):
        print_short_help(parser)
        parser.exit()



if __name__=='__main__':
    shelp="""Usage: showtime.py [options] ncfile
           [--logger-level LOGGER_LEVEL] [--logger-file LOGGER_FILE]

Show time axis dates of netcdf file.
Show time axis dates of netcdf file.

Options:
  -h, --help            show this help message and exit
  -t TNAME, --time=TNAME
                        name of the time axis variable
  -v VNAME, --var=VNAME
                        name of the variable from which to guess time
  -s TSLICE, --slice=TSLICE
                        time slice (examples: "2", ":-2:4")
  -m, --minmax          show min/max only (default: False)
  -n NCOL, --ncol=NCOL  number of comlumns for output (default: 5)
  -f FORMAT, --format=FORMAT
                        date format (default: %Y-%m-%d %H:%M:%S)
"""
#    shelp="""Options:
#  -h, --help            show a reduced help
#  --long-help           show an extended help
#  --cfgfile=CFGFILE     Configuration file [default: "config.cfg"]
#
#  Global configuration options:
#     Important general configuration options
#
#    --mars=MARS         complete configuration of mars. [default:
#                        'MANGA-V8.11']
#    --mars-version=MARS_VERSION
#                        MARS version. [default: 'V8.11']
#"""
    print opt2rst(shelp)


