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

from __future__ import absolute_import
import fnmatch, math, os, re, shutil

__author__ = 'Jonathan Wilkins'
__email__ = 'wilkins@actimar.fr'
__doc__ = '''
File utilities
==============

This module provides various file related features:
    - filesystem traversal with depth support
    - file search, wildcard or regex based
    - file rollover (backup)
    - size parsing and formatting
    - directory creation without error on existing directory
'''


def mkdirs(d):
    '''
    Create a directory, including parents.

    :Params:
        - **d**: the directory, or list of directories, that may be created

    :Return:
        - **created**: For a single directory: d if directory has been created, '' otherwise
          (already exists). For a list of directories, the list of directories which have been created.
    '''
    if not isinstance(d, basestring):
        return [dd for dd in d if mkdirs(dd)]
    if d and not os.path.exists(d):
        os.makedirs(d)
        return d
    return ''


def mkfdirs(f):
    '''
    Create a file directory, including parents. This may be used before writing to a file
    to ensure the parent directories exists.

    :Params:
        - **f**: the file, or list of files, for which the directory may be created

    :Return:
        - **created**: For a single file: f directory if it has been created, '' otherwise
          (already exists). For a list of files, the list of f directories which have been created
          were created
    '''
    if not isinstance(f, basestring):
        return [os.path.dirname(ff) for ff in f if mkfdirs(ff)]
    return mkdirs(os.path.dirname(f))


def rollover(filepath, count=1, suffix='.%d', keep=True, verbose=False):
    '''
    Make a rollover of the specified file. Keep a certain number of backups
    of a file by renaming them with a suffix number.

    :Params:
        - **filepath**: the file to make a backup of
        - **count**: maximum number of backup files
        - **suffix**: suffix to use when renaming files, must contain a '%d' marker which will be used to mark backup number
        - **keep**: whether to keep existing file in addition to the backup one

    :Return: True if a backup occured, False otherwise (count is 0 or filepath does not exists)

    '''
    if not count > 0: return False
    if not os.path.exists(filepath): return False
    fnt = '%s%s'%(filepath, suffix)
    for i in range(count - 1, 0, -1):
        sfn = fnt%(i)
        dfn = fnt%(i + 1)
        if os.path.exists(sfn):
            if os.path.exists(dfn):
                os.remove(dfn)
                if verbose: print 'rollover remove %s'%(dfn)
            os.rename(sfn, dfn)
            if verbose: print 'rollover rename %s -> %s'%(sfn, dfn)
    dfn = fnt%(1)
    if os.path.exists(dfn):
        os.remove(dfn)
        if verbose: print 'rollover remove %s'%(dfn)
    if keep:
        shutil.copy(filepath, dfn)
        if verbose: print 'rollover copy %s -> %s'%(filepath, dfn)
    else:
        os.rename(filepath, dfn)
        if verbose: print 'rollover rename %s -> %s'%(filepath, dfn)
    return True


_sort_size_dict = lambda sd: sorted(sd.items(), lambda a, b: cmp(a[1],b[1]))

# Binary units : 1 kibioctet (Kio) = 2^10 = 1024
_size_units = {
    'K':2**10, 'M':2**20, 'G':2**30,
    'T':2**40, 'P':2**50, 'E':2**60,
    'Z':2**70, 'Y':2**80,
}
_sorted_size_units = _sort_size_dict(_size_units)

# SI units : 1 kilooctet (Ko) = 10^3 = 1000
_si_size_units = {
    'K':10**3, 'M':10**6, 'G':10**9,
    'T':10**12, 'P':10**15, 'E':10**18,
    'Z':10**21, 'Y':10**24,
}
_sorted_si_size_units = _sort_size_dict(_si_size_units)

_strfsize_doc_sorted_units = ', '.join(map(lambda s:s[0], _sorted_size_units))

def strfsize(size, fmt=None, unit=None, si=False, suffix=True):
    '''
    Format a size in bytes using the appropriate unit multiplicator (Ko, Mo, Kio, Mio)

    :Params:

        * **size**:
            the size in bytes
        * **fmt**:
            the format to use, will receive size and unit arguments, if None
            formats "%%(size).3f %%(unit)s" or "%%(size)d %%(unit)s" will be automatically used.
        * **unit**:
            use an auto determinated unit if None, or the given one among %s
        * **si**:
            whether to use SI (International System) units (10^3, ...) or binary units (2^10, ...)

    :Return: a string
    '''

    units_dict = _si_size_units if si else _size_units
    units = reversed(_sorted_si_size_units if si else _sorted_size_units)
    unit_suffix = 'o' if si else 'io'
    size = float(size)

    fmt_unit, fmt_ratio = '', 1
    if unit is None:
        for unit, threshold in units:
            if size >= threshold:
                fmt_unit, fmt_ratio = unit, threshold
                break
    else:
        unit = unit.upper().strip()
        if unit not in units_dict:
            raise ValueError('Invalid unit, must be one of: %s'%(_strfsize_doc_sorted_units))
        fmt_unit, fmt_ratio = unit, units_dict[unit]

    fmt_size = size / fmt_ratio
    if fmt is None:
        fmt = '%(size).3f %(unit)s' if float(fmt_size) % 1 else '%(size)d %(unit)s'
    if suffix:
        fmt_unit += unit_suffix
    return fmt%{'size':fmt_size, 'unit':fmt_unit}

strfsize.__doc__ %= _strfsize_doc_sorted_units

_strpsizerex = re.compile(r'(?P<number>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s*(?P<unit>%s)?(?P<usfx>io|o)?'%('|'.join(_size_units.keys())), re.IGNORECASE)

def strpsize(size, si=False):
    """Parse a size in Ko, Mo, Kio, Mio, ...

    :Params:

        - **size**: the size string (eg. "1Ko", "1Kio", "2 Mo", " 10 Go"
        - **si**: when unit does not ends with 'io' force interpretation as
                  International System units (10^3, ...) instead of binary units (2^10, ...)

    :Return: the float number of bytes
    """
    if not isinstance(size, basestring): size = '%s'%(size)
    m = _strpsizerex.match(size)
    if m:
        d = m.groupdict()
        n = float(d['number'])
        u = (d.get('unit') or '').upper()
        s = (d.get('usfx') or '').lower()
        if u:
            if s == 'io':
                r = n * _size_units[u]
            elif si:
                r = n * _si_size_units[u]
            else:
                r = n * _si_size_units[u]
        else:
            r = n
        return int(math.ceil(r))
    raise ValueError('Cannot parse size: %s'%(size))


def walk(top, topdown=True, onerror=None, followlinks=False, depth=None, onfile=None, ondir=None, _depth=0):
    '''
    New implementation of os.walk **with depth support to avoid unnecessary large scans**.
    This yield a supplementary depth value for each walk (top, dirs, nondirs, depth)

    :Params:
        - **depth**: Limit the depth of walk:
            - None: no limit
            - 0: limited to top directory entries
            - 1: limited to first directory under the top directory
            - N: limited to Nth directory under the top directory

    .. warning::
        **Do not use the _depth attribute** as it is used to track the current depth in the yield processing

    :See: :func:`os.walk` for more details on other parameters.
    '''
    if depth is not None and depth >= 0 and _depth > depth: return
    try: names = os.listdir(top)
    except os.error, err:
        if onerror is not None: onerror(err)
        return
    dirs, nondirs = [], []
    for name in names:
        path = os.path.join(top, name)
        if os.path.isdir(path):
            dirs.append(name)
            if ondir: ondir(path)
        else:
            nondirs.append(name)
            if onfile: onfile(path)
    if topdown:
        yield top, dirs, nondirs, _depth
    for name in dirs:
        path = os.path.join(top, name)
        if followlinks or not os.path.islink(path):
            for x in walk(path, topdown, onerror, followlinks, depth, onfile, ondir, _depth+1):
                yield x
    if not topdown:
        yield top, dirs, nondirs, _depth


def xfind(pattern, path=None, depth=0, files=True, dirs=False, matchall=False, abspath=True, exclude=None,
         followlinks=False, expandpath=True, onerror=None, onfile=None, ondir=None, onmatch=None):
    '''
    Find paths matching the pattern wildcard.

    :Params:
        - **pattern**: pattern or list of patterns using special characters \\*,?,[seq],[!seq] (see standard module fnmatch)
        - **path**: if not None, entries are searched from this location, otherwise current directory is used
        - **depth**: if not None, it designate the recursion limit (0 based, None for no limit, see walk function)
        - **files**: if False, file entries will not be returned
        - **dirs**: if False, directory entries will not be returned
        - **matchall**: if False, only file/directory names are evaluated, entire path otherwise
        - **abspath**: if True, returned paths are absolute
        - **exclude**: if not None, it designate a pattern or list of patterns which will be used to exclude files or directories
        - **followlinks**: if True, symbolic links will be walked (see walk function)
        - **expandpath**: if True, environment variables and special character ~ will be expanded in the passed search path

    :Example:

    >>> find('*.nc', '/path/to/data')
    ['/path/to/data/data_2010-01-01.nc', '/path/to/data/data_2010-01-02.nc', ...]

    >>> find(('*.nc', '*.grb'), '/path/to/data', depth=1, exclude=('*-01.nc', '*02.grb'))
    ['/path/to/data/data_2010-01-02.nc', '/path/to/data/grib/data_2010-01-01.grb', ...]

    '''
    if not isinstance(pattern, (list, tuple)): pattern = (pattern,)
    if not isinstance(exclude, (list, tuple)): exclude = (exclude,) if exclude is not None else tuple()
    if not path: path = '.'
    if expandpath: path = os.path.expanduser(os.path.expandvars(path))
    if path.endswith(os.path.sep): path = path[:-1]
    for r, d, f, n in walk(path, topdown=True, followlinks=followlinks, depth=depth, onerror=onerror, onfile=onfile, ondir=ondir):
        e = []
        if files: e.extend(f)
        if dirs: e.extend(d)
        for n in e:
            s = matchall and os.path.join(r, n) or n
            if any((fnmatch.fnmatch(s, p) for p in pattern)):
                if any((fnmatch.fnmatch(s, x) for x in exclude)): continue
                f = os.path.join(r, n)
                if abspath: f = os.path.abspath(f)
                if onmatch: onmatch(f)
                yield f


def find(*args, **kwargs):
    '''Build a list from the :func:`xfind` generator'''
    return list(xfind(*args, **kwargs))


def xefind(regex, path=None, depth=0, files=True, dirs=False, matchall=False, abspath=True, exclude=None,
          followlinks=False, expandpath=True, onerror=None, onfile=None, ondir=None, onmatch=None,
          getmatch=False, rexflags=None, xrexflags=None):
    '''
    Find paths matching the regex regular expression.

    :Params:
        - **regex**: the file regular expression
        - **path**: if not None, entries are searched from this location, otherwise current directory is used
        - **depth**: if not None, it designate the recursion limit (0 based, None for no limit, see walk function)
        - **files**: if False, file entries will not be returned
        - **dirs**: if False, directory entries will not be returned
        - **matchall**: if False, only file/directory names are evaluated, entire path otherwise
        - **abspath**: if True, returned paths are absolute
        - **exclude**: if not None, it designate a regular expression which will be used to exclude files or directories
        - **getmatch**: if True, return a list of (path, match_object) couples
        - **followlinks**: if True, symbolic links will be walked (see walk function)
        - **regexflags**: if not None, it will be used as regex compile flags
        - **xregexflags**: if not None, it will be used as exclude regex compile flags
        - **expandpath**: if True, environment variables and special character ~ will be expanded in the passed search path

    :Example:

        >>> find('.*\.nc', '/path/to/data')
        ['/path/to/data/data_2010-01-01.nc', '/path/to/data/data_2010-01-02.nc', ...]

        >>> filelist = find('data_([0-9]{4})-([0-9]{1,2})-([0-9]{1,2})\.nc', 'data', getmatch=True, abspath=False)
        >>> for filepath, matchobj in filelist:
        >>>     print filepath, ':', matchobj.groups()
        data/data_2010-01-1.nc : ('2010', '01', '1')
        data/data_2010-01-10.nc : ('2010', '01', '10')
    '''
    if not path: path = '.'
    if expandpath: path = os.path.expanduser(os.path.expandvars(path))
    if path.endswith(os.path.sep): path = path[:-1]
    if rexflags is not None: x = re.compile(regex, rexflags)
    else: x = re.compile(regex)
    if exclude and xrexflags is not None: X = re.compile(exclude, xrexflags)
    elif exclude: X = re.compile(exclude)
    else: X = None
    for r, d, f, n in walk(path, topdown=True, followlinks=followlinks, depth=depth, onerror=onerror, onfile=onfile, ondir=ondir):
        e = []
        if files: e.extend(f)
        if dirs: e.extend(d)
        for n in e:
            s = matchall and os.path.join(r, n) or n
            m = x.match(s)
            if m:
                if X and X.match(s): continue
                f = os.path.join(r, n)
                if abspath: f = os.path.abspath(f)
                if onmatch: onmatch(f)
                yield getmatch and (f,m) or f


def efind(*args, **kwargs):
    '''Build a list from the :func:`xefind` generator'''
    return list(xefind(*args, **kwargs))


def tfind(regex, path=None, fmt='%Y-%m-%dT%H:%M:%SZ', min=None, max=None, group=None, getdate=False, getmatch=False, xmin=False, xmax=True, **kwargs):
    '''
    Find timestamped paths (e.g. files having a date string in their paths)

    :See: func:`find` for **regex**, **path** and **kwargs** arguments.

    The regex regular expression must define at least one group which describe the date string location in paths.

    :Params:

        - **fmt**: (python) date format
        - **min**: minimum date filter: a datetime object or a date string in fmt format. None means no max date filtering.
        - **max**: maximum date filter: a datetime object or a date string in fmt format. None means no max date filtering.
        - **group**: the regex group(s) number(s) or name(s): one or a list of string or integer. None means all groups.
        - **xmin**: if True, min is exclusive
        - **xmax**: if True, max is exclusive

    The group(s) can be specified either by their number or name. These group will be concatenated to
    form the date that will be parsed.

    :Examples:

    Assuming we are lokking for the follwing files:

        - path/to/data/data_2010-01-01T00H.nc
        - path/to/data/data_2010-01-01T12H.nc
        - path/to/data/data_2010-01-02T00H.nc
        - path/to/data/data_2010-01-02T12H.nc

    The commands below will have the same result:

        >>> items = tfind('data_(.*)\.nc', 'path/to', '%Y-%m-%dT%HZ', depth=2)
        >>> items = tfind('data_(....-..-..T..Z)\.nc', 'path/to/data', '%Y-%m-%dT%HZ')

    Same but more precise / advanced examples:

        >>> items = tfind('data_([0-9]{4}-[0-9]{4}-[0-9]{4}T[0-9]{2}Z)\.nc', 'path/to/data', '%Y%m%dT%HH')
        >>> items = tfind('data_([0-9]{4})-([0-9]{4})-([0-9]{4})T([0-9]{2})Z\.nc', 'path/to/data', '%Y%m%d%H')
        >>> items = tfind('(data)_(?P<y>[0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2})Z\.nc', 'path/to/data', '%Y%m%d%H', group=('y',3,4,5)))

    :Return:

        Depending on getdate and getmatch, a list in the form:

        - If getdate=False and getmatch=False: [path1, path2, ...]
        - If getdate=False and getmatch=True:  [(path1, match1), (path2, match1), ...]
        - If getdate=True  and getmatch=False: [(path1, datetime1), (path2, datetime2), ...]
        - If getdate=True  and getmatch=True:  [(path1, matchobj1, datetime1), (path2, matchobj2, datetime2), ...]

    '''
    # If min or max are strings, parse them to datetimes
    if isinstance(min, basestring): min = datetime.datetime.strptime(min, fmt)
    if isinstance(max, basestring): max = datetime.datetime.strptime(max, fmt)
    # Find paths, getting matched groups
    items = efind(regex, path, getmatch=True, **kwargs)
    # If group is not specified (or 0, meaning all groups), add concatenation of all groups
    if not group:
        items = list((i[0], i[1], ''.join(i[1].groups())) for i in items)
    # Else add concatenation of named/indexed groups
    else:
        if isinstance(group, (basestring, int)): group = (group,)
        m2s = lambda m: ''.join(isinstance(g, basestring) and m.groupdict()[g] or m.group(g) for g in group)
        items = list((i[0], i[1], m2s(i[1])) for i in items)
    # Convert matched groups to parsed datetime
    items = list((i[0], i[1], datetime.datetime.strptime(i[2], fmt)) for i in items)
    # Filter by min/max datetimes
    if min:
        if xmin: items = filter(lambda i: i[2] > min, items)
        else: items = filter(lambda i: i[2] >= min, items)
    if max:
        if xmax: items = filter(lambda i: i[2] < max, items)
        else: items = filter(lambda i: i[2] <= max, items)
    # Sort by dates
    items = sorted(items, lambda a,b: cmp(a[2], b[2]))
    # Remove non-requested fields
    if not getdate and not getmatch:
        items = list(i[0] for i in items)
    else:
        get = [0]
        if getmatch: get.append(1)
        if getdate: get.append(2)
        if not getdate or not getmatch:
            items = list(tuple(i[g] for g in get) for i in items)
    return items


