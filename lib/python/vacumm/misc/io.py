# -*- coding: utf8 -*-
"""In/Output tools"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2017)
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
import copy
import logging.handlers
import os, gc, glob, logging, re, sys
from collections import OrderedDict
from traceback import format_exc
from warnings import warn
import warnings

import numpy as N, MV2, cdms2, cdtime
import pylab as P, operator
from _geoslib import Point, LineString, Polygon
from matplotlib.collections import LineCollection, PolyCollection
from mpl_toolkits.basemap import Basemap

from ..__init__ import VACUMMError

__all__ = ['list_forecast_files', 'NcIterBestEstimate', 'NcIterBestEstimateError', 'NcFileObj',
    'ncfind_var', 'ncfind_axis', 'ncfind_obj', 'ncget_var', 'ncread_var', 'ncread_files', 'ncread_best_estimate',
    'ncget_grid', 'ncget_time','ncget_lon','ncget_lat','ncget_level', 'ncmatch_obj',
    'ncget_axis', 'netcdf3', 'netcdf4', 'ncread_axis', 'ncread_obj', 'ncget_fgrid',
    'grib_read_files', 'nccache_get_time', 'grib2nc', 'grib_get_names',
    'Shapes', 'XYZ', 'XYZMerger', 'write_snx', 'ColoredFormatter', 'Logger', 'TermColors'
    'NcIterTimeSlice']
__all__.sort()

MA = N.ma
MV = MV2


class ColPrinter(object):
    """
    Class to print formatted columns with header and frame

    :Params:

        - **columns**: A list of column descriptions like [['Year',6,'%i'],...]
          where the fist element is the title
          of the column, the second is the width of the
          column, and the last is the format of data.
        - **file**, optional: Output to file instead of stdout.
        - **align**, optional: Text alignment of titles in header, in ('left','right','center').
        - **left,right,top,bottom**, optional: String to draw part of a frame
        - **frame**: Shortcut to left=right='|',  and top=bottom='-'
        - **headsep**, optional: Header separator. If empty or None, not printed.
        - **colsep**, optional: Column separator.
        - **print_header**, optional: Printer once object created.

    .. Note::

        Obviously, you must use monospaced fonts

    :Example:

        >>> p = ColPrinter([['Year',6,'%i'],['Value',7,'%4.1f']], headsep='=',align='center',frame=True)
        ------------------
        |  Year   Value  |
        |================|
        >>> p(2000,25.9)
        | 2000   25.9    |
        >>> p(2001,23.7)
        | 2001   23.7    |
        >>> p.close()
        ------------------
    """
    def __init__(self,columns,headsep=True,align='left',
        left=None,right=None,top=None,bottom=None,frame=None,colsep=' ',
        print_header=True,file=None):
        # Inputs
        if not isinstance(columns,(list,tuple)):
            raise TypeError, 'Input must be a list or tuple of 3-element lists'
        elif len(columns) == 3 and type(columns[0]) is type('s')  and \
            type(columns[1]) is type(1) and type(columns[2]) is type('s'):
            columns = [columns,]
        if align not in ['left','right','center']:
            align = 'left'
        if frame is True:
            left = right = '|'
            top = bottom = '-'
        for bb in 'left','right','top','bottom':
            if eval(bb) is not None:
                setattr(self,'_'+bb,eval(bb))
            else:
                setattr(self,'_'+bb,'')
        self._colsep = str(colsep)
        if self._colsep == '':
            self._colsep = ' '
        if self._colsep != ' ':
            self._colsep = ' '+self._colsep+' '
        if headsep in ['',None,False]:
            self._headsep = None
        elif headsep is True:
            self._headsep = '-'
        else:
            self._headsep = str(headsep)

         # Loop on columns
        headers = []
        self._width = []
        self._fmt = []
        self._ncol = len(columns)
        for col in columns:
            if not isinstance(col,(list,tuple)) or \
              len(col) != 3 or type(col[0]) is not type('s') or type(col[1]) != type(1) or \
              type(col[2]) is not type('s') or col[2].find('%') < 0:
                raise TypeError, 'This description of column must contain the column title (string), the column width (int) and the string format for data (string): %s' % col
            width = col[1]
            if width < len(col[0]):
                width = len(col[0])
            self._width.append(width)
            if align == 'left':
                headers.append(col[0].ljust(width))
            elif align == 'right':
                headers.append(col[0].rjust(width))
            else:
                headers.append(col[0].center(width))
            self._fmt.append(col[2])

        # File case
        if file is not None:
            self._fid = open(file,'w')
        else:
            self._fid = None

        # Print header
        self._header = self._colsep.join(headers)
        if print_header:
            self.header()
        self._bottom__printed = False

    def __print(self,line):
        if isinstance(self._fid,file):
            self._fid.write(line+'\n')
        else:
            print line

    def close(self):
        """Print bottom and close the file (if file mode) """
        self.bottom()
        if isinstance(self._fid,file) and not self._fid.closed:
            self._fid.close()

    def __left(self):
        if len(self._left):
            return self._left+' '
        return ''

    def __right(self):
        if len(self._right):
            return ' '+self._right
        return ''

    def top(self):
        """ Print top frame """
        if len(self._top):
            n = len(self._header)+len(self.__left())+len(self.__right())
            self.__print((self._top * n)[:n])

    def header(self):
        """ Print header (stop + header text + header separator) """
        self.top()
        self.__print(self.__left()+self._header+self.__right())
        self.headsep()

    def headsep(self):
        """ Print header separator line """
        if self._headsep is not None:
            n = len(self._header) + \
              len(self.__left()) - len(self._left) + \
              len(self.__right()) - len(self._right)
            self.__print(self._left+(self._headsep * n)[:n]+self._right)

    def bottom(self):
        """ Print bottom frame """
        if self._bottom__printed: return
        if len(self._bottom):
            n = len(self._header)+len(self.__left())+len(self.__right())
            self.__print((self._bottom * n)[:n])
        self._bottom__printed = True


    def __call__(self,*data):
        """Print data in columns

        Arguments refer to the data.
        The number of arguments must be the same as the number of columns.
        """
        if len(data) != self._ncol:
            raise TypeError, 'You must give a number of data equal to the number of columns (%i): %i' \
            % (self._ncol,len(data))

        dds = []
        for i,dd in enumerate(data):
            dds.append((self._fmt[i] % dd).ljust(self._width[i]))
        self.__print(self.__left()+self._colsep.join(dds)+self.__right())

col_printer = ColPrinter # compat



def list_forecast_files(filepattern, time=None, check=True,
    nopat=False, patfreq=None, patfmtfunc=None, patmargin=None, verbose=False, sort=True):
    """Get a list of forecast files according to a file pattern

    :Params:

        - **filepattern**: It can be either:

            - a global matching pattern (``"file??.nc"``),
            - a date pattern (``"file%Y-%m-%d.nc"``),
            - an url (``"http://site.net/file.nc"``),
            - a list of files.

        - **time**: A time selector (``('2000', '2001', 'co')``).

          .. warning::
              This argument is *mandatory* if ``filepattern`` is a date pattern,
              and *not used* if ``filepattern`` is of another type.

        - **check**, optional: Check if local files exist.
        - **nopat**, optional: Never consider that input patterns have date patterns.
        - **patfreq**, optional: Frequency of files to generate file names for each date
           when ``filepattern`` is a date pattern.
        - **patfmtfunc**, optional: Function to use in place of
           :func:`~vacumm.misc.atime.strftime` to generate file names.
           It must take as arguments a date pattern and a CDAT component time.
        - **sort**, optional: If True, files are sorted alphabetically after being listed;
          if a callable function, they are sorted using this function (``files=sort(files)``).

          .. warning:: Files are sorted alphabetically by default!

    :Examples:

        >>> 'Prefered way'
        >>> list_forecast_files('mrsPRVMR_r0_%Y-%m-%d_00.nc', ('2010-08-06', '2010-08-15'))
        >>> list_forecast_files('http://www.ifremer.fr/data/mrsPRVMR_r0_%Y-%m-%d_00.nc', ('2010-08-06', '2010-08-15'))
        >>> list_forecast_files('mrsPRVMR_r0_%Y-%m-%d_*.nc', ('2010-08-06', '2010-08-15'))

        >>> 'Possible way'
        >>> list_forecast_files('mrsPRVMR_r0_2010-05-??_00.nc')
        >>> list_forecast_files(['mrsPRVMR_r0_2010-05-??_00.nc', 'mrsPRVMR_r0_2010-05-??_60.nc'])

        >>> 'Just ot filter in existing files'
        >>> list_forecast_files(['mrsPRVMR_r0_2010-05-06_00.nc', 'mrsPRVMR_r0_2010-05-07_00.nc'])

        >>> 'Simple conversion to list'
        >>> list_forecast_files('http://www.ifremer.fr/data/mrsPRVMR_r0_2010-05-06_00.nc')
    """
    sfp = str(filepattern)
    if len(sfp)>300: sfp = sfp[:300]+'...'
    if verbose:
        print 'Guessing file list using:'
        print '   filepattern: %s'%sfp
        print '   time selector: %s'%(time, )

    # A list of file
    if isinstance(filepattern, list):
        files = []
        for filepat in filepattern:
            files.extend(list_forecast_files(filepat, time=time, check=check,
                patfreq=patfreq, patfmtfunc=patfmtfunc, patmargin=patmargin,
                verbose=False, nopat=nopat, sort=False))


    # Date pattern
    elif not nopat and has_time_pattern(filepattern):

        from atime import pat2freq,IterDates, strftime, is_interval, pat2glob, filter_time_selector
        if isinstance(time, cdms2.selectors.Selector):
            seltime = filter_time_selector(time, noslice=True)
            if seltime.components():
                _, comps = split_selector(seltime) # FIXME: include positional components
                for i, comp in comps:
                    itv = comp.spec
                    if not is_interval(itv) and is_time(itv):
                        itv = (itv, itv, 'ccb')
                    if i==0:
                        time = itv
                    else:
                        time = itv_union(itv, time)
            else:
                time = None

        if not is_interval(time):
            if is_time(time):
                time = (time, time, 'ccb')
            else:
                raise ValueError('Your file pattern contains date pattern (like "%Y"), '
                    'so you must provide a valid absolute time interval such as (date1,date2,"co")'
                    ' or at least a valid single date')
        time = tuple(time)
        patfmtfunc = patfmtfunc if callable(patfmtfunc) else strftime

        # Guess the minimal frequency
        lmargin = 1
        if patfreq is None:
            patfreq = pat2freq(filepattern)
            if verbose: print 'Detected frequency for looping on possible dates: '+patfreq.upper()
        if not isinstance(patfreq, tuple):

            #  Reform
            patfreq = (1, patfreq)

            # Guess left margin when possible
            gfiles = glob.glob(pat2glob(filepattern))
            gfiles.sort()
            if gfiles<2:
                lmargin = 1
            elif not glob.has_magic(filepattern):
                date0 = date1 = None
                tu = 'seconds since 200'
                for i in xrange(len(gfiles)-1):
                    try:
                        date0 = strptime(gfiles[i], filepattern).torel(tu)
                        date1 = strptime(gfiles[i+1], filepattern).torel(tu)
                    except:
                        continue
                    if (date0.value>=reltime(time[0], tu).value
                        or date1.value<=reltime(time[1], tu).value): break
                if None not in [date0, date1]:
                    dt = adatetime(date1)-adatetime(date0)
                    if dt.seconds!=0:
                        lmargin = comptime('2000').add(dt.seconds, cdtime.Seconds).torel(
                            patfreq[1]+' since 2000').value
                    else:
                        lmargin = comptime('2000').add(dt.days, cdtime.Days).torel(
                            patfreq[1]+' since 2000').value
                else:
                    lmargin = 1
             #FIXME: make it work with date+glob patterns
            else:
                lmargin = patfreq[0]


#        # Add margin to time interval
#        if patmargin is None:
#
#            # Guess margin from a second time pattern (like %Y_%%Y.nc)
#            fp2 = patfmtfunc(filepattern, now())
#            if has_time_pattern(fp2):
#                for fn in glob(pat2glob(fp2):
#
#
#        elif not isinstance(patmargin, tuple):
#            patmargin = (1, patmargin)

        # Make a loop on dates
        itertime = (round_date(time[0], patfreq[1], 'floor'), time[1])
        itertime = add_margin(itertime, (lmargin-1, patfreq[1]), False)
        iterdates = IterDates(itertime, patfreq,
            closed = len(time)==3 and time[2][1]=='c' or True)
        files = []
        for date in iterdates:
            file = patfmtfunc(filepattern, date)
            if '://' in file:
                files.append(file)
            elif check or glob.has_magic(file):
                files.extend(glob.glob(file))
            else:
                files.append(file)

    # Simple remote file or file object
    elif isinstance(filepattern, cdms2.dataset.CdmsFile) or '://' in filepattern:

        files = [filepattern]

    # Glob pattern
    else:
        if check or glob.has_magic(filepattern):
            files = glob.glob(filepattern)
        else:
            files = [filepattern]

    # Unique
    files = list(set(files))

    # Count
    if verbose:
        if not files:
            print 'No file found with this file pattern "%s" and time interval %s'%(filepattern, time)
        else:
            print 'Found %i files'%len(files)

    # Sort files
    if sort:
        key = lambda x: getattr(x, 'id', x)
        if callable(sort):
            files = sort(files, key=key)
        else:
            files.sort(key=key)

    return files


def ncfind_var(f, id, ignorecase=True, regexp=False, **kwargs):
    '''
    Find a variable in a netcdf file using :func:`ncfind_obj`

    '''
    nfo = NcFileObj(f)
    f = nfo.f
    res =  ncfind_obj(f, id, ignorecase=ignorecase, regexp=regexp,
                      ids=f.listvariables(), **kwargs)
    del nfo
    return res

def ncfind_axis(f, specs, ignorecase=True, regexp=False, **kwargs):
    '''
    Find an axis in a netcdf file using :func:`ncfind_obj`

    '''
    nfo = NcFileObj(f)
    f = nfo.f
    res =  ncfind_obj(f, specs, ignorecase=ignorecase, regexp=regexp,
                      ids=f.listdimension(), **kwargs)
    del nfo
    return res

def ncfind_obj(f, specs, ignorecase=True, regexp=False, ids=None, searchmode=None, **kwargs):
    '''
    Find a variable or an axis in netcdf file using a name, list of names
    or matching attributes such as standard_name, long_name and units.

    Objects are checked using :func:`ncmatch_obj`.
    It first checks the standard_name, then the names (ids), the axis, and finally
    the long_names and units.

    :Example:

        >>> f = cdms2.open('temp.nc')
        >>> ncfind_obj(f, 'temp')
        >>> ncfind_obj(f, ['temperature','temp'])
        >>> ncfind_obj(f, ('temperature','TEMP'), ignorecase=False)
        >>> ncfind_obj(f, dict(standard_name="sea_surface_temperature"))
        >>> ncfind_obj(f, 'lon')

    :Params:

        - **f**: A cdms2.dataset.CdmsFile.
        - **name**: A string or list of string to look for,
          or a dictionary with keys "name", "standard_name", "long_name", 'units' and 'axis'.
        - **ignorecase**, optional: Ignore name case when searching variable.
        - **regexp**, optional: Interpret long_names and units as regular expressions
          that must be compiled.
        - **searchmode**, optional: Search order when ``specs`` is a dictionary
          and not a OrderedDict. It defaults
          to ``None`` or ``'snlua'`` which means:
          *standard_name -> name -> long_name -> units -> axis* (first letters).
          If ``name`` is an OrderedDict, it simply acts as a filter to restrict search.

    :Return:

        The first matching object name, or None if not found.

    '''
    nfo = NcFileObj(f)
    f = nfo.f

    # Searched terms
    withdict = isinstance(specs, dict)
    standard_name = long_name = units = axis = None
    def check_list(refname):
        if refname in specs and not is_iterable(specs[refname]):
            specs[refname] = [specs[refname]]
    if withdict: # Using a dictionary

        specs = specs.copy()

        # Useful functions
        def check_aliases(refname, *aliases):
            for alias in aliases:
                if alias in specs:
                    if refname not in specs:
                        specs[refname] = specs[alias]
                    del specs[alias]
        re_flag = re.I if ignorecase else 0
        def check_regexp(refname):
            if refname in specs and regexp: # Regular expression
                specs[refname] = [re.compile(x, re_flag).search
                                  for x in specs[refname] if x is not None]

        # Id
        check_aliases('id', 'ids', 'name', 'names')
        check_list('id')

        # Standard name
        check_aliases('standard_name', 'standard_names')
        check_list('standard_name')

        # Axis
        if 'axis' in specs:
            specs['axis'] = specs['axis'].upper()

        # Long_name
        check_aliases('long_name', 'long_names')
        check_list('long_name')
        check_regexp('long_name')

        # Units
        check_list('units')
        check_regexp('units')

    else: # Using a list of ids
        specs = {"id": specs}
        check_list('id')

    # Order and filter the search
    specs = filter_search(specs, searchmode)

    # Loop on targets to match attributes
    if ids is None:
        ids = f.listvariables()+f.listdimension()
    elif isinstance(ids, basestring):
        ids = [ids]
    for att, val in specs.items(): # loop on attribute types
        for id in ids: # Loop on targets
            if match_atts(f[id], {att:val}, id=True, ignorecase=ignorecase):
                break
        else:
            continue
        break
    else: # Not found
        id = None

    nfo.close()

    return id

def filter_search(specs, searchmode=None):
    """Order and filter the search"""
    # - get order
    all_keys = ['id', 'standard_name', 'long_name', 'units', 'axis']
    all_keys0 = [key[0] for key in all_keys]
    if searchmode is None:
        if isinstance(specs, OrderedDict):
            searchmode = ''.join([key[0] for key in specs.keys()])
        else:
            searchmode = ''.join(all_keys0)
    searchmode = searchmode.replace('n', 'i')
    keys = []
    for key0 in searchmode:
        key = all_keys[all_keys0.index(key0)]
        if key0 in all_keys0 and key in specs:
            keys.append(key)
    # - reorder specs
    return OrderedDict([(key, specs[key]) for key in keys])

def ncmatch_obj(obj, id=None, standard_name=None,
        long_name=None, units=None, axis=None, ignorecase=True,
        searchmode=None, **kwargs):
    """Check if an MV2 object (typicaly from a netcdf file) matches names, standard_names, etc


    It first checks the standard_name, then the names (ids), the axis, and finally
    the long_names and units.

    :Params:

        - **obj**: A MV2 array.
        - **standard_name**, optional: List of possible standard_names.
        - **id**, optional: Name (id) of this array, wich defaults to the id attribute.
        - **axis**, optional: Axis type, as one of 'x, 'y', 'z', 't'.
        - **long_name**, optional: List of possible long_names or callable expression
          (such as regular expression method).
        - **units**, optional: Same as ``long_names`` but for units.

    :Example:

        >>> ncmatch_obj(sst, standard_name='sea_surface_temperature', id=['sst'])
        >>> import re
        >>> ncmatch_obj(sst, long_name=re.compile('sea surface temp').match)
    """
    # Format
    search = OrderedDict()
    if id is None and 'name' in kwargs:
        id = kwargs['name']
    for key in ('standard_name', 'id', 'axis', 'long_name', 'units'):
        val = locals()[key]
        if val is None:
            continue
        search[key] = val
        if key=='axis':
            search[key] = val if not isinstance(key, list) else val[0]
            continue
        search[key] = val

    # Order and filter the search
    search = filter_search(search, searchmode)

    # Check attributes
    return match_atts(obj, search, ignorecase)


def ncget_var(f, *args, **kwargs):
    '''
    Get a variable object as returned by :meth:`cdms2.dataset.CdmsFile.getVariable`
    which is equivalent to ``f[varname]``.

    :Return: A cdms2.fvariable.FileVariable or None if not found.

    :See: :func:`ncfind_var()`
    '''
    nfo = NcFileObj(f)
    f = nfo.f
    oldvname = args[0]
    vname = ncfind_var(f, *args, **kwargs)
    if vname is None:
        raise IOError('Variable not found %s in file %s'%(oldvname, f.id))
    var =  f.getVariable(vname)
    del nfo
    return var


def ncread_obj(f, name, *args, **kwargs):
    """Read an arbitrary netcdf object (axis or variable)"""
    # Inits
    nfo = NcFileObj(f)
    f = nfo.f
    ignorecase = kwargs.get('ignorecase', True)
    ncname = ncfind_obj(f, name, ignorecase=ignorecase)
    if ncname is None:
        del nfo
        raise IOError('Object not found %s in file %s'%(name, f.id))
    if ncname in f.variables:
        obj = ncread_var(f, ncname, *args, **kwargs)
    else:
        obj = ncread_axis(f, ncname, *args, **kwargs)
    del nfo
    return obj

def ncread_axis(f, name, select=None, ignorecase=True, mode='raise'):
    """Read a 1D axis

    .. note:: Please use :func:`ncread_var` to read 2D axes.

    :Params:

        - **mode**, optional: if ``'raise'`` raises an :exc:`IOError`
          if not found, else returns ``None``.
    """
    # Inits
    nfo = NcFileObj(f)
    f = nfo.f
    ncname = ncfind_axis(f, name, ignorecase=ignorecase)
    if ncname is None:
        del nfo
        if mode=='raise':
            raise IOError('Axis not found %s in file %s'%(name, f.id))
        return
    axis = f.getAxis(ncname)
    atts = get_atts(axis)
    axis = axis.clone()
    set_atts(axis, atts)
    del nfo
    return axis


def ncread_var(f, vname, *args, **kwargs):
    """Read a variable in a netcdf file and some more

    In addition to a simple ``f(vname, *args, **kwargs)```:

        - ``vname`` can be a list of var names, and it takes the first
          one found, ignoring the case by default.
        - If a variable is on a grid that is stored as curvilinear
          but is rectangular in real, it convert its grid to a rectanguar grid

    If variabe is not found, it raises

    :Params:

        - **f**: File descriptor.
        - **vname**: Variable name(s) (see :func:`ncfind_var`).
        - **ignorecase**, optional: Case insensitive search for the name of variable.
        - Other arguments and keywords are passed to ``f``.
        - **atts**: Dictionary of attributes to apply.

        - **squeeze**, optional: A single argument (or a list of them) interpreted
          as a squeeze specification passed to :func:`~vacumm.misc.misc.squeeze_variable`,
          to squeeze out singleton axes.
        - **torect**, optional: If True, try to convert output grid to rectanguar
          using :func:`~vacumm.misc.grid.misc.curv2rect`.
        - **mode**, optional: if ``'raise'`` raises an :exc:`IOError`
          if not found, else returns ``None``.

    :Example:

        >>> var = ncread_var(f, ['u', 'u2'], lon=(45, 47), atts={'id':'U'})

    """
    # Inits
    nfo = NcFileObj(f)
    f = nfo.f
    ignorecase = kwargs.pop('ignorecase', True)
    torect = kwargs.pop('torect', True)
    grid = kwargs.pop('grid', None)
    kwgrid = kwfilter(kwargs, 'grid')
    samp = kwargs.pop('samp', None)
    atts = kwargs.pop('atts', None)
    mode = kwargs.pop('mode', 'raise')
    searchmode = kwargs.pop('searchmode', None)
    squeeze = kwargs.pop('squeeze', False)
    if not isinstance(squeeze, list):
        squeeze = [squeeze]
    oldvname = vname

    # Find
    vname = ncfind_var(f, vname, ignorecase=ignorecase)
    if vname is None:
        del nfo
        if mode=='raise':
            raise VACUMMError('Variable not found %s in file %s'%(oldvname, f.id))
        return

    # Read
    var = f(vname, *args, **kwargs)

    # Extra stuff
    var = _process_var(var, torect, samp, grid, kwgrid, squeeze, atts)

    del nfo
    return var

def _process_var(var, torect, samp, grid, kwgrid, squeeze, atts):
    '''
    - **samp**, optional: Undersample rate as a list of the same size as
      the rank of the variable. Set values to 0, 1 for no undersampling.
    - **torect**, optional: If True, try to convert output grid to rectanguar
      using :func:`~vacumm.misc.grid.misc.curv2rect`.
    - **grid**, optional: A grid to regrid the variable on.
    - **grid_<keyword>**, optional: ``keyword`` is passed to
      :func:`~vacumm.misc.grid.regridding.regrid`.
    - **squeeze**, optional: Argument passed to :func:`squeeze_variable` to squeeze out singleton axes.
            - Extra kwargs are used to refine the **selector** initialized with ``select``.
    - **atts**: attributes dict (or list of attributes dict for each varname)
    '''
    # To rectangular grid?
    if torect:
        curv2rect(var, mode="none")

    # Undersample
    if samp is not None:
        if len(samp)!=var.ndim:
            del nfo
            raise VACUMMError('Sampling keyword ("samp") must have'
                ' a size equal to the rank of the variable (%i)'%var.ndim)
        for i, ss in enumerate(samp):
            if samp==0 or samp==1:
                samp[i] = slice(None)
            elif not isinstance(ss, slice):
                samp[i] = slice(None, None, ss)
        var = var(*samp)

    # Regrid
    if grid is not None:
        var = regrid2d(var, grid, **kwgrid)

    # Squeeze
    if squeeze:
        for ss in squeeze:
            if ss is False: break
            var = squeeze_variable(var, ss)
            if ss in [True, None]:
                break

    # Attributes
    if atts is not None:
        set_atts(var, atts)

    return var

def ncget_axis(f, checker, ids=None, ro=False, checkaxis=False, **kwargs):
    """Get an axis in a netcdf file by searching all axes and variables

    If ``checker`` is a list, dict or tuple, :func:`ncfind_axis`
    is called directly to search for the axis within the file.

    :Param:

        - **checker**: Can be either

            - A generic name such as 'x' or 'lon',
            - A function to check that an object is an axis.
              of appropriate type (such as :func:`~vacumm.misc.axes.islon`).
              This function must accept the 'ro' keyword ('readonly').
            - An argument to :func:`ncfind_axis`: list, dict, tuple.

        - **ids**, optional: A list of ids to focus search.

    :Return: The axis or None if not found
    """

    nfo = NcFileObj(f)
    f = nfo.f
    if isinstance(checker, basestring):
        checker = get_checker(checker)
    elif isinstance(checker, (list, tuple, dict)):
        axid = ncfind_obj(f, checker, ids=ids, checkaxis=checkaxis, ro=ro, **kwargs)
        if axid is None: return
        axis = f[axid]
        atts = get_atts(axis)  # to fix a cdms bug
        axis = axis.clone()
        set_atts(axis, atts)
        nfo.close()
        return axis

    # Targets
    dimids = f.listdimension()
    varids = f.listvariables()
    if ids is None:
        ids = varids+dimids
    elif isinstance(ids, basestring):
        ids = [ids]

    # Loop on targets
    axis = None
    for id in ids:
        if checker(f[id], checkaxis=checkaxis, ro=True):
            axis = f[id]
            break
#            if id in varids:
#                return f(id)
#            return f.getAxis(id)
        elif id in varids:
            for i, aid in enumerate(f[id].listdimnames()):
                if checker(f[aid], checkaxis=checkaxis, ro=True):
                    axis = f[id].getAxis(i)
                    break
            else:
                continue
            break
    if axis is None:
        return
    if hasattr(axis, 'clone'):
        atts = get_atts(axis)  # to fix a cdms bug
        axis = axis.clone()
        set_atts(axis, atts)
    else:
        axis = axis()
    del nfo
    checker(axis, checkaxis=checkaxis, ro=ro)
    return axis


def ncget_lon(f, ids=None, checkaxis=False, ro=False):
    """Get longitude axis of a netcdf file

    :Params:

        - **f**: Netcdf file name or object.
        - **ids**, optional: List of ids to help searching.
    """
    return ncget_axis(f, islon, ids, checkaxis=checkaxis, ro=ro)

def ncget_lat(f, ids=None, checkaxis=False, ro=False):
    """Get latitude axis of a netcdf file

    :Params:

        - **f**: Netcdf file name or object.
        - **ids**, optional: List of ids to help searching.
    """
    return ncget_axis(f, islat, ids, checkaxis=checkaxis, ro=ro)

def ncget_time(f, ids=None, checkaxis=False, ro=False):
    """Get time axis of a netcdf file

    :Params:

        - **f**: Netcdf file name or object.
        - **ids**, optional: List of ids to help searching.
    """
    return ncget_axis(f, istime, ids, checkaxis=checkaxis, ro=ro)

def ncget_level(f, ids=None, checkaxis=False, ro=False):
    """Get level axis of a netcdf file

    :Params:

        - **f**: Netcdf file name or object.
        - **ids**, optional: List of ids to help searching.
    """
    return ncget_axis(f, islevel, ids, checkaxis=checkaxis, ro=ro)

def ncget_grid(f, ids=None, torect=False):
    """Get a grid of a netcdf file

    :Params:

        - **f**: Netcdf file name or object.
        - **ids**, optional: List of ids to help searching.
    """
    nfo = NcFileObj(f)
    f = nfo.f

    if ids is None:
        ids = f.listvariables()
    elif isinstance(ids, basestring):
        ids = [ids]
    grid = None
    for id in ids:
        grid = f[id].getGrid()
        if grid is not None :
#            grid = grid.clone()
            break
    del nfo
    if torect:
        grid = curv2rect(grid, mode="none")
    return grid

def ncget_fgrid(f, gg):
    """Get the file grid that matches a transient grid or variable

    Matching is checked using ids of longitudes and latitudes.

    :Params:

        - **f**: file name or object.
        - **gg**: cdms2 grid or variable with a grid.

    :Return: A :class:`FileGrid` instance or ``None``
    """
    if f is None or gg is None: return
    if isgrid(f): return f

    # From a variable
    if cdms2.isVariable(gg):
        gg = gg.getGrid()
        if gg is None: return

    # Current grid
    if isinstance(gg, tuple):
        lon, lat = gg
    else:
        lon = gg.getLongitude()
        lat = gg.getLatitude()
    if lon is None or lat is None: return
    lonid = getattr(lon, '_oldid', lon.id)
    latid = getattr(lat, '_oldid', lat.id)

    # Loop on file grids
    from vacumm.misc.io import NcFileObj
    nfo = NcFileObj(f)
    for fgrid in nfo.f.grids.values():
        flon = fgrid.getLongitude()
        flat = fgrid.getLatitude()
        flonid = getattr(flon, '_oldid', flon.id)
        flatid = getattr(flat, '_oldid', flat.id)
        if (lonid, latid)==(flonid, flatid):
            nfo.close()
            return fgrid
    nfo.close()


_nccache_time = {}
def nccache_get_time(f, timeid=None, ro=False):
    """Get a time axis from cache or netcdf file

    A time axis not in cache is read using :func:`ncget_time`,
    them stored in cache.

    :Params:

        - **f**: File object or name.
        - **timeid**, optional: Single or list of time ids for :func:`ncget_time`.

    Example:

        >>> taxis = nccache_get_time('myfile.nc', ['time','time_counter'])
    """
    # Check cache
    # - file name
    fname = f if isinstance(f, basestring) else f.id
    fname = os.path.realpath(fname)
    # - from cache
    if fname in _nccache_time:
        return _nccache_time[fname]

    # Read it
    taxis = ncget_time(f, ids=timeid, ro=ro, checkaxis=True)
    _nccache_time[fname] = taxis
    return taxis

class NcIterTimeSlice(object):
    """Basic netcdf file iterator with a fixed slice"""
    def __init__(self, files, tslice=None, timeid=None, keepopen=False, autoclose=True):
        self.i = 0
        if isinstance(files, basestring): files = [files]
        self.nfiles = len(files)
        self.files = files
        self.nfonext = None
        if tslice is None:
            tslice = slice(None)
        self.tslice = tslice
        self.autoclose = autoclose
        self.keepopen = keepopen
        self.tslices = []
        self.timeid = timeid

    def __iter__(self):
        return self

    def next(self, verbose=False):

        # Last iteration
        if self.i == self.nfiles:
            self.close()
            raise StopIteration

        # Open current file
        if self.nfonext is not None: # from next one
            f = self.nfonext.f
            if not self.keepopen: self.nfo.close()
            self.nfo = self.nfonext
            self.nfonext = None
        else: # first time used
            self.nfo = NcFileObj(self.files[self.i])
            f = self.nfo.f

        # Get time axis
        if hasattr(f, '_vacumm_timeid'):
            self.timeid = f._vacumm_timeid
        taxis = nccache_get_time(f, timeid=self.timeid, ro=True)
        if taxis is None:
            return f, None
        self.timeid = taxis.id
        if verbose:
            print taxis

        # Real slice
        tslice = slice(*self.tslice.indices(len(taxis)))

        # Finalize
        self.i += 1
        return f, self.tslice

    def close(self):
        """Close file descriptors that can be closed"""
        if not self.keepopen:
            if self.nfonext:
                self.nfonext.close()
            if self.nfo:
                self.nfo.close()

class NcIterBestEstimate(object):
    """Iterator on netcdf forecast files

    This class is useful for reading the best estimate of netcdf forcast files.

    :Params:

        - **files**: A list of netcdf files.
        - **toffset**, optional: An integer or tuple of (<num>, '<units>')
          to skip the first time steps of each files.
        - **timeid**, optional: Time id. If ``None``, it is guessed
          using :meth:`guess_timeid`.
        - **tslices**, optional: A list of time slices
          (typically taken from a previous loop on file),
          to prevent guessing them.
        - **keepopen**, optional: Keep all file descriptor open,
          else close those who can be closed once no longer used.
        - **autoclose**: Deprecated.

    :Iterator: At each iteration, it returns ``f,tslice``

        - ``f``: the file descriptor (may be closed),
        - ``tslice``: the time slice

            - A :class:`slice` instance.
            - ``None`` if not time is found (thus no slice to perform).
            - ``False``: if nothing to read at all.

    :Example:

    >>> for f, tslice in NcIterBestEstimate(ncfiles, toffset=(1,'day')):
    ...     if tslice is False or time is None: continue
    ...     var = f(time=tslice)
    """
    def __init__(self, files, time=None, toffset=None, timeid=None, tslices=None, keepopen=False, autoclose=True, id=None):
        self.i = 0
        if isinstance(files, basestring): files = [files]
        self.nfiles = len(files)
        self.files = files
        self.nfonext = None
        self.seltime = time #[time] if isinstance(time,list) else time
        self.tslices = [] if tslices is None else tslices
        if toffset is None: toffset = 0
        self.toffset = toffset
        self.timeid = timeid
        self.autoclose = autoclose
        self.keepopen = keepopen
        from atime import add
        self.add = add
        if id is None:
            id = str(self.toffset)+str(self.seltime)+str(files)
            try:
                import md5
                id = md5.md5(self.id).digest()
            except:
                pass
        self.id = id

    def __iter__(self):
        return self

    def next(self, verbose=False):

        # Last iteration
        if self.i == self.nfiles:
            self.close()
            raise StopIteration

        # Check cache of time slices
        if len(self.tslices)>self.i:
            self.nfo = NcFileObj(self.files[self.i])
            f = self.nfo.f
            tslice = self.tslices[self.i]
            self.i += 1
            return f, tslice

        # Open current file
        if self.nfonext is not None: # from next one
            f = self.nfonext.f
            if not self.keepopen: self.nfo.close()
            self.nfo = self.nfonext
            self.nfonext = None
        else: # first time used
            self.nfo = NcFileObj(self.files[self.i])
            f = self.nfo.f

        # Check file cache
        if not hasattr(f, '_vacumm_nibe_tslices'):
            f._vacumm_nibe_tslices = {}

        # Get time axis
        if hasattr(f, '_vacumm_timeid'):
            self.timeid = f._vacumm_timeid
        taxis = nccache_get_time(f, timeid=self.timeid, ro=True)
        if taxis is None:
            return f, None
        self.timeid = taxis.id
        if verbose:
            print taxis

        # Get base time slice of current file
        ijk = tsel2slice(taxis, self.seltime, asind=True, nonone=True)
        if ijk is False:
            return self.empty() # nothing
        i, j, k = ijk

        # Offset
        ctimes = taxis.asComponentTime()
        if isinstance(self.toffset, tuple):
            subseltime = (self.add(ctimes[0], *self.toffset), ctimes[-1], 'cc')
            subtaxis = taxis.subAxis(*ijk)
            ijo = subtaxis.mapIntervalExt(subseltime)
            if ijo is None or ijo[2]==-1: return self.empty() # nothing
            i += ijo[0]
        else:
            i = max(i, self.toffset)
            if i>=j: return self.empty() # nothing
            subtaxis = None

        # Truncate to next file
        if self.i+1 != self.nfiles:

            # Start of slice
            if subtaxis is None:
                subtaxis = taxis.subAxis(*ijk)
                ct0 = subtaxis.subAxis(0, 1).asComponentTime()[0]
            else:
                ct0 = ctimes[i]

            # End of slice
            self.nfonext = NcFileObj(self.files[self.i+1])
            taxisnext = nccache_get_time(self.nfonext.f, timeid=self.timeid, ro=True)
            if isinstance(self.toffset, tuple):
                ct1 = taxisnext.subAxis(0, 1).asComponentTime()[0]
                ct1 = self.add(ct1, *self.toffset)
                bb = 'co'
            else:
                if self.toffset>=len(taxisnext): # Next file too short for offset
                    ct1 = ctimes[-1]
                    bb ='cc'
                else: # First step starting from offset
                    ct1 = taxisnext.subAxis(self.toffset, self.toffset+1).asComponentTime()[0]
                    bb = 'co'
            subseltime = (ct0, ct1, bb)
            ijo = subtaxis.mapIntervalExt(subseltime)
            if ijo is None or ijo[2]==-1: return self.empty() # nothing
            io, jo, ko = ijo
            j = min(j, i+jo)
        del ctimes

        # Finalize
        self.i += 1
        tslice = slice(i, j)
        self.tslices.append(tslice)
        f._vacumm_nibe_tslices[self.id] = tslice
        return f, tslice

    def empty(self):
        """Nothing to read from this file"""
        self.tslices.append(None)
        self.i += 1
        self.nfo.f._vacumm_nibe_tslices[self.id] = None
        return self.nfo.f, False

    def close(self):
        """Close file descriptors that can be closed"""
        if not self.keepopen:
            if self.nfonext:
                self.nfonext.close()
            if self.nfo:
                self.nfo.close()


class NcIterBestEstimateError(VACUMMError):
    pass


class NcFileObj(object):
    """Simple class to properly manage file object or name

    :Examples:

        >>> nfo = NcFileObj('myfile.nc')
        >>> nfo.f('sst')
        >>> nfo.close() # or del nfo: close file descriptor

        >>> f = cdms2.open(path)
        >>> nfo = NcFileObj(f)
        >>> nfo.f('sst')
        >>> del nfo # or nfo.close(): here has no effect (descriptor still open)
        >>> f.close()
        >>> nfo = NcFileObj(f)
        >>> nfo.close() # close file descriptor

    """
    def __init__(self, ncfile, mode='r'):
        if isinstance(ncfile, NcFileObj):
            self.type = ncfile.type
            self.f = ncfile.f
        elif isinstance(ncfile, basestring):
            self.type = 'path'
            self.f = cdms2.open(ncfile, mode)
        elif hasattr(ncfile, '_status_'):
            self.type = ncfile._status_
            if self.type == 'closed':
                self.f = cdms2.open(ncfile.id, mode)
            else:
                self.f = ncfile
        else:
            raise IOError('Unknown file type %s (not a file name or a netcdf file object)'%ncfile)
    def isclosed(self):
        return self.type == 'closed'
    def ispath(self):
        return self.type == 'path'
    def isopen(self):
        return self.type == 'open'
    def close(self):
        if self.type in ['closed', 'path'] and self.f._status_ == 'open':
            self.f.close()
    __del__ = close

def ncread_best_estimate(filepattern, varname, *args, **kwargs):
    """Read the best estimate of a variable through a set of netcdf forecast files

    .. warning:: This function is deprecated.
        Please use :func:`ncread_files` switching the first
        two argument.

    This is equivalent to::

        ncread_files(varname, filepattern, *args, **kwargs)
    """
    return ncread_files(varname, filepattern, *args, **kwargs)

def ncread_files(filepattern, varname, time=None, timeid=None, toffset=None, select=None,
    atts=None, samp=None, grid=None, verbose=False, ignorecase=True, torect=True,
    squeeze=False, searchmode=None, nibeid=None, sort=True, nopat=False, patfreq=None,
    patfmtfunc=None, check=True, bestestimate=True, **kwargs):
    """Read the best estimate of a variable through a set of netcdf files

    .. warning:: Files are listed using function :func:`list_forecast_files`.
        Please read its documentation before using current function.

    :Examples:

        >>> var = ncread_files("r0_2010-%m-%d_00.nc", 'xe',
            ('2010-08-10', '2010-08-15', 'cc'), samp=[2, 1, 3])
        >>> var = ncread_files("http://www.net/r0_2010-%m-%d_00.nc", 'xe',
            ('2010-08-10', '2010-08-15', 'cc'),
            timeid='TIME', toffset=(1, 'day'))
        >>> var = ncread_files("r0_2010-??-??_00.nc", 'xe',
            select=dict(lon=(-10,-5), z=slice(23,24)), grid=smallgrid)
        >>> xe, sst = ncread_files("myfiles*.nc", [('xe', 'sla'),('sst','temp'),'u'])

    :Params:

        - **varname**: Name of the netcdf variable to read.

            - If a simple name, it reads this variable.
            - If a list of names, it reads them all.
            - If a list of list of names, each variable is searched for
              using the sublist of names.

        - **filepattern**: File pattern. See :func:`list_forecast_files`
          for more information.
        - **time**, optional: Time selector. This keyword is *mandatory*
          if ``filepattern`` has date patterns.
        - **toffset**: Skip the first time steps. See :class:`NcIterBestEstimate`
          for more information.
        - **select**, optional: An additional selector for reading
          the variable. It can be a dictionary or a :class:`~cdms2.selectors.Selector`
          instance (see :func:`~vacumm.misc.misc.create_selector`).
        - **atts**: attributes dict (or list of attributes dict for each varname)
          (see :func:`ncread_var`.)
        - **samp**, optional: Undersample rate as a list of the same size as
          the rank of the variable. Set values to 0, 1 for no undersampling.
        - **grid**, optional: A grid to regrid the variable on.
        - **grid_<keyword>**, optional: ``keyword`` is passed to
          :func:`~vacumm.misc.grid.regridding.regrid`.
        - **timeid**, optional: Time id (otherwise it is guessed).
        - **ignorecase**, optional: Ignore variable name case (see :func:`ncfind_var`).
        - **torect**, optional: If True, try to convert output grid to rectanguar
          using :func:`~vacumm.misc.grid.misc.curv2rect` (see :func:`ncread_var`).
        - Extra kwargs are used to refine the **selector** initialized with ``select``.
        - **squeeze**, optional: Argument passed to :func:`ncread_var`
          to squeeze out singleton axes.
        - **searchmode**, optional: Search order (see :func:`ncfind_obj`).
        - **sort/nopat/patfreq/patfmtfunc/check**, optional: These arguments are passed to
          :func:`list_forecast_files`.

    :Raise: :class:`NcIterBestEstimateError` in case of error.
    """
    # Get the list of files
    ncfiles = list_forecast_files(filepattern, time, sort=sort, nopat=nopat,
        patfreq=patfreq, patfmtfunc=patfmtfunc, check=check)
    if len(ncfiles)==0:
        raise NcIterBestEstimateError('No valid file found with pattern: %s'%filepattern)
    single = not isinstance(varname, list)
    varnames = [varname] if single else varname
    if verbose:
        print 'Reading best estimate variable(s): ', ', '.join([str(v) for v in varnames]), '; time:', time
        print 'Using files:'
        print '\n'.join([getattr(fn, 'id', fn) for fn in ncfiles])

    # Some inits
    nvar = len(varnames)
    if isinstance(atts, dict):
        atts = [atts]
    atts = broadcast(atts, nvar)
    allvars = [[] for iv in xrange(nvar)]
    kwgrid = kwfilter(kwargs, 'grid')
    # - base selector
    selects = broadcast(select, nvar)
    selectors = [create_selector(s, **kwargs) for s in selects]
    # - iterator on files
    if not toffset and not bestestimate and (isinstance(time, slice) or time is None):
        iterator = NcIterTimeSlice(ncfiles, time, timeid=timeid)
    else:
        iterator = NcIterBestEstimate(ncfiles, time, timeid=timeid, toffset=toffset, id=nibeid)
    # - undepsampling
    if samp is not None:
        samp = [0 for s in samp if s==0 or not isinstance(s, int)]
        samp = [slice(None, None, s) for s in samp]
    # - output grid
    if grid is not None:
        from grid.regridding import regrid2d
    # - time
    time_units = None
    newgrid = None
    tvars = [False]*len(varnames) # vars with time?
    itaxes = {}

    # Loop on files
    for ifile, (f, tslice) in enumerate(iterator):

        # Refine selector specs with time slice
        kwseltime = {iterator.timeid:tslice} if iterator.timeid is not None and \
            isinstance(tslice, slice) and not tslice==slice(None) else {}
#        seltime = selector(**kwseltime)
        taxis = None

        # Loop on variables
        for iv, vn in enumerate(varnames):

            # Refine this selector
            seltime = selectors[iv](**kwseltime)

            # Find variable name
            oldvn = vn
            vn = ncfind_var(f, vn, ignorecase=ignorecase)
            if vn is None:
                if verbose: print 'Skipping file %s for %s variable not found'%(f.id, oldvn)
                continue

            # Check time
            if f[vn] is None:
                continue
            withtime = iterator.timeid is not None and iterator.timeid in f[vn].getAxisIds()
            if withtime:
                itaxes[iv] = f[vn].getOrder().find('t')
                tvars[iv] = True
                if not tslice:
                    if verbose: print 'Skipping file %s for %s variable because time slice not compatible'%(f.id, oldvn)
                    continue
                sel = seltime # with time
            else:
                sel = selectors[iv] # no time

            # Infos
            if verbose:
                print 'Processing file no', ifile, ' ', f, ', variable:', vn, ', time slice :', tslice
                if withtime:
                    if taxis is None: taxis = f[vn].getTime()
                    ctimes = taxis.asComponentTime()
                    print '  Available:', ctimes[0], ctimes[-1]
                    del ctimes

            # Read the variable
            if verbose:
                print '  Selecting:', sel
            try:
                var = ncread_var(f, vn, sel, ignorecase=True, torect=torect, squeeze=squeeze,
                                 grid=grid, samp=samp, searchmode=searchmode,
                                 atts=atts[iv] if atts is not None and iv<len(atts) else None)
                if verbose:
                    print '  Loaded:', var.shape
            except Exception, e:
                if verbose: print 'Error when reading. Skipping. Message: \n'+format_exc()#e.message
                continue


            # Fix time units (that may vary between files)
            if iterator.timeid is not None and withtime:
                this_time_units = f[iterator.timeid].units
                if time_units is None:
                    time_units = this_time_units
                elif not are_same_units(this_time_units, time_units):
                    try:
                        ch_units(var, time_units)
                    except:
                        continue

            # Update
            if withtime or ifile==0: # append first time or variables with time
                allvars[iv].append(var)
            if True not in tvars and ifile==1: # no time for all variables
                break

        gc.collect()

        # Read only one file if no variable with time
        if ifile==0 and True not in tvars:
            break

    iterator.close()

    # Concatenate
    from misc import MV2_concatenate
    for iv in xrange(nvar):

        # Check
        if len(allvars[iv])==0:
            raise VACUMMError('No valid data found using varname(s): %s, '
                'filepattern: %s, time: %s'%(varnames[iv], filepattern, time))

        # Reorder and merge
        allvars[iv] = MV2_concatenate(allvars[iv], axis=itaxes.get(iv, 0), copy=False)

    return allvars[0] if single else allvars



def grib_get_names(gribfile):
    '''
    Return a list of a grib file parameter unique names (using grib message's shortName).
    '''
    import pygrib
    names = []
    with pygrib.open(gribfile) as g:
        for i in xrange(g.messages):
            m = g.read(1)[0]
            if m.shortName not in names:
                names.append(m.shortName)
    return names


def grib_read_files(
    filepattern, varname, time=None, select=None,
    torect=None, samp=None, grid=None, squeeze=None, atts=None,
    verbose=False, **kwargs):
    """
    Read cdms2 variables through one or a set of grib files.

    :Examples:

        >>> vardict = grib_read_files("r0_2010-%m-%d_00.grb", 'u',
                ('2010-08-10', '2010-08-15', 'cc'), samp=[2, 1, 3])
        >>> vardict = grib_read_files("r0_2010-??-??_00.grb", dict(shortName:'u'),
                select=dict(lon=(-10.0,-5.0), lat=slice(100,200)), grid=smallgrid)
        >>> vardict = grib_read_files("myfiles*.grb", [dict(shortName=['u', 'u10']), dict(shortName=['v','v10'])])

    :Params:

        - **filepattern**: must be:
            - File pattern. See :func:`list_forecast_files` for more information.
            - One or more string(s) of the files(s) to be processed. string(s) may contain wildcard characters.
        - **varname**: Name of the grib variable(s) to read.
            - If a simple name, it reads this variable **using the grib message's shortName**.
            - If a list of names, it reads them all.
            If a name is a dict, then it is used as grib selector in which case
            the user should not specify selectors which may interfer with the select keyword
            (see :func:`~pygrib.open.select`).
        - **time**, optional: Time selector for files and data. This keyword is *mandatory*
          if ``filepattern`` has date patterns.

        - **select**, optional: An additional selector applied after data have been loaded.
          It can be a dictionary or a :class:`~cdms2.selectors.Selector`
          instance (see :func:`~vacumm.misc.misc.create_selector`).

        - **torect**, optional: If True, try to convert output grid to rectanguar
          using :func:`~vacumm.misc.grid.misc.curv2rect` (see :func:`ncread_var`).
        - **samp**, optional: Undersample rate as a list of the same size as
          the rank of the variable. Set values to 0, 1 for no undersampling.
        - **grid**, optional: A grid to regrid the variable on.
        - **grid_<keyword>**, optional: ``keyword`` is passed to
          :func:`~vacumm.misc.grid.regridding.regrid`.
        - **squeeze**, optional: Argument passed to :func:`ncread_var` to squeeze out singleton axes.
        - **atts**: attributes dict (or list of attributes dict for each varname)

        - *verbose*: function to be called for logging (sys.stderr if True,
            disabled with False)

    :Return:

        If varname is a list of names or dicts:
        - a dict of loaded variables as :class:`cdms2.tvariable.TransientVariable`
          this dict keys are are filled with the corresponding varname value if it is a string, or wiht
          the loaded message's shortName/name/parameterName.
        Else:
        - the loaded variable as :class:`cdms2.tvariable.TransientVariable`

    """
    import datetime, glob, numpy, os, sys, time as _time, traceback
    import cdms2, pygrib
    from axes import create_lat, create_lon
    from atime import create_time, datetime as adatetime
    from grid import create_grid, set_grid
    if not verbose:
        verbose = lambda s:s
    if verbose and not callable(verbose):
        verbose = lambda s: sys.stderr.write(('%s\n')%s)
    # List of variables
    single = not isinstance(varname, (list,tuple))
    varnames = [varname] if single else varname
    verbose(
        'grib_read_files:\n'
        '  filepattern: %s\n'
        '  time: %s\n'
        '  varname: %s'
    %(filepattern, time, '\n- '.join(['%r'%(v) for v in varnames])))
    # List of files
    if isinstance(filepattern, basestring):
        files = list_forecast_files(filepattern, time)
    else:
        if not isinstance(filepattern, (list, tuple)):
            filepattern = (filepattern,)
        files = tuple(f for l in map(lambda p: glob.glob(p), filepattern) for f in l)
    if len(files)==0:
        raise Exception('No valid file found with pattern %r and time %r'%(filepattern, time))
    verbose('number of matching files: %s'%(len(files)))
    #verbose('- %s'%('\n- '.join(files)))
    if time:
        time = map(adatetime, time[:2])
    vardict = dict()
    # Load grib data
    for f in files:
        verbose('file: %s'%(f))
        with pygrib.open(f) as g:
            for n in varnames:
                kw = n if isinstance(n, dict) else dict(shortName=n)
                st = _time.time()
                ms = g.select(**kw)
                verbose('  select: %s message%s matching preselection %r (took %s)'%(
                    len(ms), 's' if len(ms)>1 else '', kw, datetime.timedelta(seconds=_time.time()-st)))
                for m in ms:
                    st = _time.time()
                    # use provided special datetime object if present
                    if m.validDate:
                        dt = m.validDate
                    # use validityDate exposed as YYYYMMDD and validityTime exposed as HHMM (or HMM or HH or H)
                    elif m.validityDate != None and m.validityTime != None:
                        dt = '%s%04d00'%(m.validityDate, m.validityTime) # pad validityTime and add 00 seconds
                    # or use dataDate & dataTime & forecastTime ??
                    else:
                        raise Exception('Don\'t know how to handle datetime for message:\n%r'%(m))
                    if isinstance(dt, basestring):
                        dt = datetime.datetime.strptime(dt, '%Y%m%d%H%M%S')
                    if time and (dt < time[0] or dt >= time[1]):
                        continue
                    if m.gridType == 'regular_ll':
                        latitudes,longitudes = m.distinctLatitudes, m.distinctLongitudes
                    else:
                        latitudes,longitudes = m.latlons()
                    kn = n
                    if isinstance(kn, dict):
                        kn = n.get('shortName', n.get('name', n.get('parameterName', None)))
                    if not kn: kn = m.shortName
                    if not kn in vardict: vardict[kn] = []
                    vardict[kn].append(dict(
                        datetime=dt,
                        latitudes=latitudes, longitudes=longitudes,
                        values=m.values
                    ))
                    verbose('    message name: %r, shortName: %r, datetime: %s, gridType: %r, latitude%s, longitude%s (took %s)'%(
                        m.name, m.shortName, dt, m.gridType, latitudes.shape, longitudes.shape, datetime.timedelta(seconds=_time.time()-st)))
                    del m
                del ms
    # Transform loaded data into cdms2 variable
    kwgrid = kwfilter(kwargs, 'grid')
    for n,p in vardict.iteritems():
        if not p:
            vardict[n] = None
            continue
        p = sorted(p, lambda a,b: cmp(a['datetime'], b['datetime']))
        time = create_time([pp['datetime'] for pp in p])
        lat = create_lat(p[0]['latitudes'])
        lon = create_lon(p[0]['longitudes'])
        var = cdms2.createVariable(
            [pp['values'] for pp in p],
            id='_'.join(n.split()), long_name=n,
        )
        var.setAxis(0, time)
        set_grid(var, create_grid(lon, lat))
        vatts = atts=atts[iv] if atts is not None and iv<len(atts) else None
        if select:
            var = var(**select)
        var = _process_var(var, torect, samp, grid, kwgrid, squeeze, atts)
        vardict[n] = var
    # Return variable or dict of variables
    return vardict.values()[0] if (single and vardict) else vardict


def grib2nc(filepattern, varname):
    '''
    ***Currently for test purpose only***
    '''
    varlist = grib_read_files(filepattern, varname, verbose=True)
    if ncoutfile:
        print>>sys.stderr, 'Writing to netcdf file:', ncoutfile
        netcdf3()
        if os.path.exists(ncoutfile):
            print>>sys.stderr, 'File already exists:', ncoutfile
            sys.exit(1)
        f = cdms2.open(ncoutfile, 'w')
        try:
            for n,v in varlist.iteritems():
                if v is None:
                    print>>sys.stderr, '  %r not found'%(n)
                    continue
                print>>sys.stderr, '  %r (%r)'%(v.id, n)
                f.write(v)
        finally: f.close()



class CachedRecord:
    """Abstract class for managing cached records

    ** It cannot be used by itself **

    The following class variables must be defined:

    - _cache_file: cache file name
    - _time_range: ('2000-12-01 12:00','2005-10-01','co') OR (1, 'day', 'cc')
    - _var_info: (('lon','Longitude', 'degrees east', None),('hs','Wave Height', 'm', (0., 20)),)
    - _station_info: (['Start Bay', '2007-12-25'],['Long Bay', '2004-05-18'])
    - _dt: (1800, cdtime.Seconds)
    - _time_units: 'seconds since 2000-01-01'

    The following methods must be defined:

    - _load_from_source_:
    """
    _verbose_level = 1
    _missing_value = 1.e20

    def __init__(self, time_range, **kwargs):

        for key, val in kwargs.items():
            key = '_'+key
            if hasattr(self, key):
                setattr(self, key, val)

        self._cache_file = self._cache_file %vars()
        self._update_mode = time_range in [None, 'update']
        if self._update_mode:
            time_range = 'all'
        self._time_range = self._get_time_range_(time_range)

        _var_names = [sn for sn, ln, un, vr in _var_info]
        _station_info = [[name, comptime(time_origin)] for name,  time_origin in _buoy_info]

        self._vars = {}
        self.check_cache()


    def _print_(self, text, level=2):
        if self._verbose_level >= level:
            if level == 1:
                text = '!WARNING! '+text
            text = '[%s] %s' %(self.buoy_type, text)
            print text

    def _warn_(self, text):
        self._print_(text, level=1)

    def show_variables(self):
        """Print available variables"""
        print 'Available variables:'.upper()
        for n,v in self._vars.items():
            print '%s: %s [%s]' % (n,v.long_name,v.units)

    def get(self,var_name, time_range=None):
        """Get a variable

        - **var_name**: Name of the variable

        Return: A 1D :mod:`MV2` variable

        See: :meth:`plot()`, :meth:`show_variables()`
        """
        # Time range
        time_range = self._get_time_range_(time_range)
        if time_range != self._time_range:
            self.load(time_range)
        assert var_name.lower() in self._vars.keys(),' Wrong name of variable ("%s"). Please use .show_variables() to list available variables'%var_name.lower()
        if time_range is None:
            args = []
        else:
            args = [time_range]
        return self._vars[var_name](*args)

    def __call__(self,*args,**kwargs):
        """Get a variable

        @see: :meth:`get()`
        """
        return self.get(*args,**kwargs)

    def plot(self,var_name,time_range=None,**kwargs):
        """Plot a variable

        - **var_name**: Name of the variable
        - *time*: Plot only within this time range (like ('2007-01-01','2007-02-01','co')
        - *show*: Show the figure [default: None]
        - All other keywords are passed to :func:`vacumm.misc.plot.curve()`

        Example:

        >>> myvar = mybuoy.plot('baro')

        See: :meth:`get()`, :meth:`show_variables()`
        """
#       # Time range
#       time_range = self._get_time_range_(time_range)
#       if time_range != self._time_range:
#           self.load(time_range)
#
        # Variable
        assert var_name.lower() in self._vars.keys(),' Wrong name of variable. Please use .show_variables() to list available variables'
        var = self._vars[var_name]
        if time_range is not None:
            var = var(time_range)

        # Keywords
        defaults = {
            'linestyle':'-',
            'marker':'.',
            'ms':4.,
            'mfc':'#00ffff',
            'title':'%s at %s buoy %s' % (var.long_name,self.buoy_type,self.buoy_id)
        }
        for att,val in defaults.items():
            kwargs.setdefault(att,val)

        # Plot
        curve(var,**kwargs)

    def save(self, file_name, mode='w', warn=True):

        if len(self._vars) == 0:
            self._warn_('No variables to save')
            return
        if file_name.endswith('.nc'):
            if not os.access(file_name, os.W_OK):
                if warn:
                    self._warn_('No write access to file:'+file_name)
                return
            f = cdms2.open(file_name, mode)
            for var in self._vars.values():
                f.write(var)
            for att in 'type', 'id', 'name':
                att = 'buoy_'+att
                setattr(f, att, getattr(self, att))
            if hasattr(self, '_url'):
                f.url = self._url
            f.close()
            if file_name == self._cache_file:
                self._print_('Cache updated')
            else:
                self._print_('Saved to '+file_name)


    def __init__(self, buoy, time_range, **kwargs):

        # Search for the buoy
        buoy_names = _channelcoast_list_(buoy)
        assert len(buoy_names), 'No buoy matching: '+buoy
        self.buoy_name = buoy_names[0]
        self.buoy_id = buoy_id = self.buoy_name.replace(' ', '')
        for buoy_name, time_origin in self._buoy_info:
            if buoy_name == self.buoy_name:
                self._time_origin = time_origin
                break
        self.buoy_type = self.__class__.__name__

        # Tuning
        for key, val in kwargs.items():
            key = '_'+key
            if hasattr(self, key):
                setattr(self, key, val)

        self._cache_file = self._cache_file %vars()
        self._update_mode = time_range in [None, 'update']
        if self._update_mode:
            time_range = 'all'
        self._time_range = self._get_time_range_(time_range)
        self._vars = {}
        self.check_cache()

    def get(self, var_name, time_range=None):
        time_range = self._get_time_range_(time_range)
        if time_range != self._time_range:
            self.load(time_range)
        return _Buoy_.get(self, var_name, time_range)

    def plot(self, var_name, time_range=None, **kwargs):
        time_range = self._get_time_range_(time_range)
        self.load(time_range)
        return _Buoy_.plot(self, var_name, time_range, **kwargs)

    def _get_time_range_(self, time_range):

        if isinstance(time_range, (list, tuple)): # Directly specified
            from .atime import time_selector
            time_range = time_selector(*time_range)

        elif time_range in  ('full', 'all', True, False): # Everything possible
            time_range = (self._time_origin, now(True))

        else :#if hasattr(self, '_time_range'): # Default specified time range
            time_range = self._time_range

        if len(time_range) == 2: # Default bound indicators
            time_range += ('co', )

        if time_range[0] < self._time_origin: # Check not too old
            time_range = (self._time_origin, ) + time_range[1:]
        return time_range

    def _check_time_range_(self, time_range, time_axis, getbounds=False):
        """Check if a time_range is included in an time axis"""

        # Format time range
        time_range = self._get_time_range_(time_range)

        # Get right time axis
        if isinstance(time_axis, dict):
            if not len(time_axis):
                res = (False, False)
                if getbounds:
                    res += (None, None)
                return res
            time_axis = time_axis.values()[0]
        if cdms2.isVariable(time_axis):
            time_axis = time_axis.getTime()

        # Get axis range
        nt = len(time_axis)
        t0, t1 = time_axis.subAxis(0, nt, nt-1).asComponentTime()

        # Convert to closed bounds
        time_range = list(time_range)
        for ibound, isign in (0, 1), (1, -1):
            if time_range[2][ibound] == 'o':
                time_range[ibound] = time_range[ibound].add(isign*self._dt[0], self._dt[1])

        # Check bounds
        res = (time_range[0] >= t0 or self._update_mode, time_range[1] <= t1)
        if getbounds:
            res += (t0, t1)
        return res

    def load(self, time_range):
        """Check if time_range is included in in-memory variables, else, check cache and load from it"""
        time_range = self._get_time_range_(time_range)
        if self._vars == {} or \
            self._check_time_range_(time_range, self._vars) != (True, True):
            self.check_cache(time_range)

    def check_cache(self, time_range=None):
        """Update the cache"""

        buoy_id = self.buoy_id

        # Time bounds
        time_range = self._get_time_range_(time_range)
        self._print_('CHECK CACHE TIME RANGE '+str(time_range))
        t0_request, t1_request, bb = time_range

        # Update or create?
        if time_range is None:raise 'a'
        if not os.path.exists(self._cache_file%vars()):  # Create
            self._print_('CREATE')
            time_range = (time_range[0], time_range[1].add(-1, cdtime.Day), 'cc')
            self._vars = self._load_from_source_(time_range[:2]+('cc', ))
            self.save(self._cache_file%vars(), warn=False)

        else: # Update

            # Check cache time range
            f = cdms2.open(self._cache_file)
            cache_time = f.getAxis('time')
            atts = get_atts(cache_time)
            cache_time = cache_time.clone()
            set_atts(cache_time, atts)
            f.close()
            t0_good, t1_good, t0_cache, t1_cache = self._check_time_range_(time_range, cache_time, getbounds=True)

            # Check first date
            vars_new = {}
            if not t0_good:
                # First date of cache is too recent
                vars_before = self._load_from_source_((t0_request, t0_cache, time_range[2][0]+'o'))
                vars_current = self._load_from_cache_('all')
                if len(vars_before):
                    for vn in self._var_names:
                        vars_new[vn] = (vars_before[vn], vars_current[vn])
            # Check last date
            if not t1_good:
                vars_after = self._load_from_source_((t1_cache, t1_request, 'o'+time_range[2][1]))
                if len(vars_after):
                    if not len(vars_new):
                        # We just append to file
                        self._print_('APPEND TO FILE')
                        self._vars = vars_after
                        self.save(self._cache_file%vars(), 'a', warn=False)
                        self._print_(' loading from cache '+str(time_range))
                        self._load_from_cache_(time_range)
                        return
                    else:
                        # We append to var before saving
                        for vn in self._var_names:
                            vars_new[vn] += (vars_after[vn], )

            # Here we completely change the cache
            if len(vars_new):

                # Merge
                self._print_('MERGING')
                for ivar, (sn, ln, un, vr) in enumerate(self._var_info):
                    self._vars[sn] = MV.concatenate(vars_new[sn])
                    self._vars[sn].id = self._vars[sn].name = sn
                    self._vars[sn].long_name = ln
                    self._vars[sn].units = un

                # Save
                self.save(self._cache_file%vars(), warn=False)

            # Update in memory variables
            if self._check_time_range_(time_range, self._vars) != (True, True):
                self._load_from_cache_(time_range)

        gc.collect()

    def _load_from_cache_(self, time_range=None):
        """Load variables from cache"""
        time_range = self._get_time_range_(time_range)
        self._vars = {}
        buoy_id  = self.buoy_id
        f = cdms2.open(self._cache_file%vars())
        for var_name in f.variables.keys():
            self._vars[var_name] = f(var_name, time_range)
        f.close()
        return self._vars




class Shapes(object):
    """A class to read shapefiles and return GEOS objects
    Inspired from basemap.readshapefile

    Here are the conversion rules from shapefile to GEOS objects :

        - Points and multipoints are interpreted as :class:`Points`.
        - Polylines are interpreted as :class:`LineString`.
        - Polygons are interpreted as :class:`Polygons`.

    :Params:

        - **input**: Refers to a shapefile or is a shapes isntance ;
          if a shapefile, it assumes that <input>.shp contains points,
          multipoints, lines or polygons, and that <input>.dbf contains their attributes.
        - **proj*, optional: A projection function to convert coordinates. It must accept
          the "inverse" keyword.
        - **m*, optional: A Basemap instance for converting for plotting.
        - *inverse*, optional: Inverset the conversion with proj .
        - **clip*, optional: If in the form ``(xmin,ymin,xmax,ymax)``,
          clips to this box ; if a polygon like argument,
          it clips to this polygon
          (see :func:`~vacumm.misc.grid.masking.polygons()` for arguments).
          If simply ``True`` and *m* is present, it clips to the bounds of *m*.
        - **min_area*, optional: Minimal area to keep a polygon
        - **samp**, optional: An integer to undersample coordinates of polygons and lines.
        - **shapetype**, optional:

            - If 0, it must only deal with points ;
            - if 1, only polylines ;
            - if 2, only polygons (conversion 1<->2 is automatic).
    """
    POINTS = POINT = 0
    LINES = LINE = 1
    POLYGONS = POLYS = POLY = POLYGON = 2
    INPUT_POINTS = 1
    INPUT_MULTIPOINTS = 8
    INPUT_POLYLINES = 3
    INPUT_POLYGONS = 5
    def __init__(self, input, m=None, proj=False, inverse=False, clip=True,
            shapetype=None, min_area=None, sort=True, reverse=True, samp=1):

        # Inits
#        if isinstance(input, list) and input: input = input[0]
        from_file = isinstance(input, str)
        if hasattr(m, 'map'): m = m.map
        default_proj = None if m is None else m
        self._m = m

        if from_file:

            # From a shapefile
            if input.endswith('.shp') and input.endswith('.dbf'):
                input = input[:-4]
            for ext in ('shp', ):#, 'dbf':
                fname = '%s.%s'%(input, ext)
                assert os.path.exists(fname), fname
            try:
                try:
                    from shapefile import Reader
                except:
                    from mpl_toolkits.basemap.shapefile import Reader
                newreader = True
                shp = Reader(input)
                input_type = shp.shapeType
            except Exception, e:
                print>>sys.stderr, 'Cannot read %s:\n%s\nTrying with shapelib'%(input, e)
                from shapelib import ShapeFile
                newreader = False
                shp = ShapeFile(input)
                input_type = shp.info()[1]
            self._prefix = input
#           dbf = dbflib.open(input)
            if default_proj and (1, 1) == default_proj(1, 1):
                default_proj = None
            self._info = []

        elif isinstance(input, (list, N.ndarray)): # From coordinates
            in_coords = input
            input_type = 5 if not len(in_coords) or in_coords[0].ndim==2 else 1
            self._info = []

        else:

            # From a Shapes (or super) instance
            in_coords = input.get_data(proj=False)
            self._m = input._m # overwrite m keyword
            default_proj = input._proj
            input_type = [Shapes.INPUT_POINTS, Shapes.INPUT_POLYLINES,
                Shapes.INPUT_POLYGONS][input._type]
            self._info = input._info

        # Get coordinates
        if from_file:
            if newreader:
                nshapes = shp.numRecords
            else:
                nshapes = shp.info()[0]
        else:
            nshapes = 1
        coords = []
        if input_type in [Shapes.INPUT_POINTS, Shapes.INPUT_MULTIPOINTS]: # A Point or MultiPoint file
            if shapetype is not None and shapetype != self.POINTS:
                raise TypeError, 'Your shape type is not point'
            self._type = self.POINTS

            # Loop on shape groups
            for iobj in xrange(nshapes):
                if from_file:
                    if newreader:
                        all_points = shp.shape(iobj).shape.points
                    else:
                        all_points = shp.read_object(iobj).vertices()
                else:
                    all_points = in_coords
                coords.extend(all_points)

            # Merge coordinates
            xy = N.asarray(coords)

#               if from_file: self._info.append(dbf.read_record(iobj))

        elif input_type in [Shapes.INPUT_POLYLINES, Shapes.INPUT_POLYGONS]: # A Polyline or Polygon file

            # Shape type
            if shapetype is not None:
                if shapetype == Shapes.POINTS:
                    raise TypeError, 'Your shape type is point, not polyline or polygon'
                else:
                    self._type = shapetype
            else:
                if input_type==Shapes.INPUT_POLYLINES:
                    self._type = Shapes.LINES
                else:
                    self._type = Shapes.POLYGONS

            # Loop on shape groups
            for iobj in xrange(nshapes):
                if from_file:
                    if newreader:
                        obj = shp.shapeRecord(iobj).shape
                        all_points = obj.points
                        if len(all_points)==0: continue
                        nparts = len(obj.parts)
                        if nparts==1:
                            all_polys = [all_points]
                        else:
                            all_polys = []
                            for ip in xrange(nparts-1):
                                all_polys.append(all_points[obj.parts[ip]:obj.parts[ip+1]])#xxxxxxxx
                    else:
                        all_polys = shp.read_object(iobj).vertices()
                else:
                    all_polys = in_coords
                coords.extend(all_polys)

            # Merge coordinates
            xy = N.concatenate(coords)

        else:
            raise TypeError, 'Input shapefile must only contains 2D shapes'

        # Bounds
        if xy.shape[0]>0:

            self.xmin = xy[:, 0].min()
            self.xmax = xy[:, 0].max()
            self.ymin = xy[:, 1].min()
            self.ymax = xy[:, 1].max()
        else:
            self.xmin = N.inf
            self.xmax = -N.inf
            self.ymin = N.inf
            self.ymax = -N.inf
        del xy
        self.xpmin = N.inf
        self.xpmax = -N.inf
        self.ypmin = N.inf
        self.ypmax = -N.inf

        # Projection
        if callable(proj):
            self._proj = proj
        elif default_proj is not None and proj is None:
            self._proj = default_proj
        elif proj is True or isinstance(proj, basestring):
#            self._proj = None
            self._proj = self._get_proj_(proj)
        else:
            self._proj = False
        self._m_projsync = None
        if callable(self._m): # same projection as of map?
            if proj is False:
                self._m_projsync = N.allclose((1, 1), self._m(1, 1))
            elif self._proj is self._m:
                self._m_projsync = True
            elif callable(self._proj) and self._proj is not self._m:
                self._m_projsync = N.allclose(self._proj(1, 1), self._m(1, 1))


        # Clipping zone with projected coordinates
        clip = self._clip_zone_(clip)

        # Convert to shapes
        self._shaper = [Point, LineString, Polygon][self._type]
        self._shapes = []
        for coord in coords:

            # Numeric array
            coord = N.asarray(coord)

            # Under sampling
            if samp > 1 and coord.shape[0] > (2*samp+1):
                coord = coord[::samp]

            # Projection
            if self._proj:
                if coord[..., 1].max()<91 and coord[..., 1].min()>-91:
                    coord[..., 1] = N.clip(coord[..., 1], -89.99, 89.99)
                coord = N.asarray(self._proj(coord[..., 0], coord[..., 1])).T

            # Convert to shape instance
            shape = self._shaper(coord)

            # Clip
            if clip:
                shapes = clip_shape(shape, clip)
            else:
                shapes = [shape]

            # Minimal area
            if min_area is not None and self._shaper is Polygon and min_area > 0.:
                shapes = filter(lambda sh: sh.area() >= min_area, shapes)

            # Store
            self._shapes.extend(shapes)


        # Final bounds
        if clip is not None or min_area:

            # Normal coordinates
            xy = self.get_xy(proj=False)
            self.xmin = xy[0].min() # min(self.xmin, xy[0].min())
            self.xmax = xy[0].max() # max(self.xmax, xy[0].max())
            self.ymin = xy[1].min() # min(self.ymin, xy[1].min())
            self.ymax = xy[1].max() # max(self.ymax, xy[1].max())
            del xy

        # Projected coordinates
        xyp = self.get_xy(proj=None)
        self.xpmin = xyp[0].min() #min(self.xpmin, xyp[0].min())
        self.xpmax = xyp[0].max() #max(self.xpmax, xyp[0].max())
        self.ypmin = xyp[1].min() #min(self.ypmin, xyp[1].min())
        self.ypmax = xyp[1].max() #max(self.ypmax, xyp[1].max())
        del xyp


        # Finalize
        if from_file and not newreader:
            shp.close()
#           dbf.close()
        # Sort by area?
        if sort:
            self.sort(reverse=reverse)
        else:
            self._sorted = 0

    def _clip_zone_(self, clip):
        """Return a projected polygon or None"""
        # From grid
        if isgrid(clip):
            lon = clip.getLongitude().getValue()
            lat = clip.getLatitude().getValue()
            clip = dict(lon=(lon.min(), lon.max()), lat=(lat.min(), lat.max()))

        # From dictionary
        if isinstance(clip, dict):
           if 'lon' in clip and 'lat' in clip:
               clip = [clip['lon'][0], clip['lat'][0], clip['lon'][1], clip['lat'][1]]

           else:
               clip = None

        # No clipping
        if clip is False: return

        # Guess or set it
        if clip is not None:

            # Normal polygon
            if clip is not True:
                return create_polygon(clip, proj=self._proj)

            # Guess from map
            if self._m is not None:
                if self._m_projsync:
                    proj = False
                    data = N.asarray([self._m.xmin, self._m.ymin, self._m.xmax, self._m.ymax])
                else:
                    xx, yy = self._m([self._m.xmin, self._m.xmax, self._m.xmax, self._m.xmin],
                        [self._m.ymin, self._m.ymin, self._m.ymax, self._m.ymax], inverse=True)
                    data = N.asarray([xx, yy])
                return create_polygon(data, proj=False)



    def clip(self, zone, copy=True, sort=True, reverse=True, **kwargs):
        """Clip to zone

        :Params:

            - **zone**: ``[xmin, ymin, xmax, ymax]``
            - **copy**, optional: If ``True``, make a copy of current instance,
              else simply rehandled the list of shapes.
            - If ``copy==True``, Other parameters are passed to the initialization
              of the new instance.
        """
        if not copy:
            zone = self._clip_zone_(zone)
            if zone is None: return self
            if self._type > 0:
                newshapes = []
                for shape in self:
                    if zone.intersects(shape):
                        intersections = shape.intersection(zone)
                        newshapes.extend(filter(lambda s: isinstance(s,  self._shaper), intersections))
                self._shapes = newshapes
                if sort: self.sort(reverse=reverse)
            return self
        return self.__class__(self, clip=zone, **kwargs)

#   def join_lines(self, poly=True):
#       # Array of end points
#       nsh = len(self)
#       ends = N.zeros((2*nsh, 3))
#       for i in xrange(nsh):
#           ends[i, 0:1] = self._shapes[i][0]
#           ends[i+1, 0:1] = self._shapes[i][-1]
#           ends[i:i+2, 2] = float(i)
#       # Distances
#       xx = ends[:, 0].reshape((nsh*2, nsh*2))
#       yy = ends[:, 1].reshape((nsh*2, nsh*2))
#       dst = (xx-xx.transpose())**2+(yy-yy.transpose())**2
#       new_shapes = []


    def __len__(self):
        return len(self._shapes)

    def __getitem__(self, key):
        return self._shapes[key]

#    def project(self, proj=True, inverse=False):
#        """Project shapes using proj
#
#        - **proj**: A Basemap instance or pure projection instance (from Basemap).
#          If True, a mercator projection is used.
#        - *inverse*: Inverse projection [default: False]
#        """
#        if proj is True:
#            proj=merc(lon=(self.xmin, self.xmax), lat=(self.ymin, self.ymax))
#        for i, shape in enumerate(self):
#            bb = shape.boundary
#            if self._type == 0:
#                bb = proj(b[0], b[1], inverse=inverse)
#            else:
#                bb[:, 0], bb[:, 1] = proj(bb[:, 0], bb[:, 1], inverse=inverse)
#            self._shapes[i] = self._shaper(bb)
#        if isinstance(proj, Basemap): proj = proj.projtran
#        self._proj.append((proj, inverse))

    def sort(self, reverse=True):
        """Sort shapes according to their surface or length

        - *reverse*: If True, greater polygons are first [default: True]
        """
        if self._type == 2:
            self._shapes.sort(cmp=lambda p0, p1: cmp(p0.area(), p1.area()), reverse=reverse)
            self._sorted = 1-2*int(reverse)
        elif self._type==1:
            self._shapes.sort(cmp=lambda p0, p1: cmp(len(p0.boundary), len(p1.boundary)), reverse=reverse)
            self._sorted = 1-2*int(reverse)
        else:
            self._sorted = 0
    def sorted(self):
        return self._sorted

    def get_type(self):
        """Return the type of shapes

        - :attr:`POINTS` = Points,
        - :attr:`LINES` = LineStrings = PolyLines
        - :attr:`POLYGONS` = Polygons
        """
        return self._type

    def is_type(self, type):
        """Check type

        :Example:

            >>> self.is_type(self.POLYS)
        """
        return self.get_type()==type

    def _get_proj_(self, proj=None):
        """Get a valid projection function to operate on shapes coordinates"""
        if proj is None: return
        if proj is True or isinstance(proj, basestring):

            if hasattr(self, '_proj') and callable(self._proj): # already projected
                return

            if proj is True and callable(self._m): # from map

                return self._m

            # using grid.basemap.get_proj
            if self.xmin>self.xmax:
                gg = None
            else:
                gg = ([self.xmin,self.xmax],
                    N.clip([self.ymin,self.ymax], -89.99, 89.99))
            kw = dict(proj=proj) if isinstance(proj, basestring) else {}
            return get_proj(gg, **kw)

        if callable(self._proj): # already projected

            if proj is False: # no projection -> project back
                return lambda x, y, inverse=False: self._proj(x, y, inverse=not inverse)

            if proj is not self._proj: #re-projection
                old_proj = proj
                def proj(x, y, inverse=False):
                    if inverse:
                        return self._proj(*(old_proj(x, y, True)+(False, )))
                    return old_proj(*self._proj(x, y, True))

            else: # no need to project
                proj = None

        return proj



    def get_shapes(self, key=None, proj=None):
        """Get the list of geos objects (polygons, etc)

        :Param:

            - **key**: A slice selector applied to the list.
            - **proj**: ``True``, or a callable to project or re-project coordinates.
        """
        # Objects to work on
        if key is None:
            shapes = self._shapes
        else:
            shapes = self._shapes[key]
        single = not isinstance(shapes, list)
        if single: shapes = [shapes]

        # Projection
        proj = self._get_proj_(proj)

        # Loop on shapes
        polys = []
        for poly in shapes:
            if proj:
                poly = proj_shape(poly, proj)
            polys.append(poly)

        if single: return polys[0]
        return polys

    def get_data(self, key=None, proj=None):
        """Get the numeric version of the list of geos objects (polygons, etc)

        :Param:

            - **key**: A slice selector applied to the list.
            - **proj**: ``True``, or a callable to project or re-project coordinates.
        """
        if not len(self): return []
        # Projection
        proj = self._get_proj_(proj)

        # Shapes to work on
        if key is None:
            shapes = self._shapes
        else:
            shapes = self._shapes[key]
        single = not isinstance(shapes, list)
        if single: shapes = [shapes]

        # Loop on shapes
        data = []
        for poly in shapes:
            xy = poly.boundary
            if proj:
                if self.is_type(self.POINTS): # Points
                    xy = N.array(proj(*xy))
                else: # Lines
                    xx, yy = proj(*xy.T)
                    xy = N.asarray([xx, yy]).T
                    del xx, yy
            data.append(xy)

        if single: return data[0]
        return data

    def get_points(self, key=None, split=True, proj=None):
        """Get all the points from all the shapes as a tuple (x,y)"""
        if not len(self): return N.array([[],[]])

        # Projection
        proj = self._get_proj_(proj)

        # Shapes to work on
        if key is None:
            shapes = self._shapes
        else:
            shapes = self._shapes[key]
        single = not isinstance(shapes, list)
        if single: shapes = [shapes]

        # Loop in shapes
        xx, yy = [], []
        for poly in shapes:
            xy = poly.boundary
            if self.is_type(self.POINTS):
                xx.append(xy[0])
                yy.append(xy[1])
            else:
                xx.extend(xy[:, 0].tolist())
                yy.extend(xy[:, 1].tolist())
        if split:
            if proj: return proj(N.asarray(xx), N.asarray(yy))
            return N.asarray(xx), N.asarray(yy)
        if proj: return N.asarray(proj(xx, yy))
        return N.asarray([xx, yy])

    def get_xy(self, key=None, proj=None):
        """Shortcut to ``get_points(split=false)``"""
        return self.get_points(split=False, key=key, proj=proj)
    xy = property(get_xy, doc="XY coordinates as a (2,npts) array")


    def resol(self, deg=True):
        """Compute the mean "resolution" of the shapes based on the first shape

        - **deg**:

            - if ``False``: return a resolution in meters has a the median distance between points
            - if ``True``: return the median distance between points as a resolution in degrees ``(xres,yres)``

        """
        if not len(self): return 0,0
        from vacumm.misc.grid.misc import resol
        x, y = self.get_xy(key=0)
        if deg and callable(self._proj): # m->deg
            dx, dy = resol((x, y), proj=False)
            x0 = x.mean()
            y0 = y.mean()
            x1, y1 = self._proj(x0+dx, y0+dx, inverse=True)
            return x1-x0, y1-y0
        elif not deg and not callable(self._proj):
            return resol((x, y), proj=True)
        return resol((x, y), proj=False)


    def get_map(self):
        """Return the associated basemap instance if set"""
        if hasattr(self._m, 'map'): return self._m.map
        return self._m

    def plot(self, select=None, ax=None, fill=None, points=False, lines=True,
            fillcolor=None, color='k', s=None, linewidth=None, m=None, show=True,
            alpha=1, autoscale=True, title=None, **kwargs):
        """Plot shapes

        :Params:

            - **select**, optional: argument for selecting shapes in the list [defaul: None].
            - **fill**, optional: Force filling (True/False), else guessed from shpe type, ie filling for polygons only [default: None]
            - **ax**, optional: Axes instance.
            - **m**, optional: :class:`~vacumm.misc.core_plot.Map` instance
              (created with :func:`~vacumm.misc.plot.map2`) or a :class:`~mpl_toolkits.basemap.Basemap` instance.
            - **points**, optional: Plots shapes as points.
            - **lines**, optional: Plot shapes as lines (if a of type :attr:`POINTS`).
            - **fill_<params>**, optional: ``<param>`` is passed to
              :class:`~matplotlib.collections.PolyCollection`.
            - **lines_<params>**, optional: ``<param>`` is passed to
              :class:`~matplotlib.collections.LineCollection` or to :class:`~matplotlib.collections.PolyCollection`.
            - **points_<params>**, optional: ``<param>`` is passed to
              :class:`~matplotlib.pyplot.scatter`.
            - **m_<params>**, optional: ``<param>`` is passed to
              :class:`~vacumm.misc.plot.map2` if ``m is True``.
            - **autoscale**, optional: Autoscale axis limits?
        """

        # Keywords
        kwpoints = kwfilter(kwargs, 'points')
        kwlines = kwfilter(kwargs, 'lines')
        kwpoints.setdefault('c', color)
        if s is not None: kwpoints.setdefault('s', s)
        if linewidth is not None:
            kwlines.setdefault('linewidth', linewidth)
        kwlines.setdefault('linestyles', 'solid')
        kwlines.setdefault('color', color)
        kwlines.setdefault('alpha', alpha)
        kwpoints.setdefault('alpha', alpha)
        kwfill = kwfilter(kwargs, 'fill')
        kwfill.update(kwlines)
        if fill is None and self.is_type(self.POLYGONS): fill = True


        # Map
        if m is True or m=='auto':
            if m!='auto' and getattr(self, '_m', None):
                m = self._m
            else:
                from core_plot import Map
                m = Map.get_current(axes=ax) or True
        if m is True:
            kwmap = kwfilter(kwargs,'m_')
            if not len(self):
                warn('No shape found, thus nothing to plot')
            else:
                if not kwmap.has_key('lon'):
                    dx = (self.xmax-self.xmin)*.05
                    kwmap['lon'] = (self.xmin-dx, self.xmax+dx)
                if not kwmap.has_key('lat'):
                    dy = (self.ymax-self.ymin)*.05
                    kwmap['lat'] = (self.ymin-dy, self.ymax+dy)
            kwmap.setdefault('res', None)
            kwmap.setdefault('proj', 'merc')
            kwmap.update(show=False, axes=ax, title=title)
            m = map2(**kwmap)
            ax = m.axes
        isbm = isinstance(m, Basemap)

        # Plot on what?
        if ax is None:
            ax = getattr(m, 'axes',  None) or P.gca()

        # Plot what?
        if points:
            xx, yy = self.get_points(split=True, proj=m)
        if not self.is_type(self.POINTS):
            data = self.get_data(select, proj=m)

        # Polygons or lines
        oo = []
        if not self.is_type(self.POINTS):
            if (fill is None and self.is_type(self.POLYGONS)) or fill is True: # Polygons
                if fillcolor is None: fillcolor=land
                for kv in dict(facecolor=fillcolor).items():
                    kwfill.setdefault(*kv)
                cc = PolyCollection(data, **kwfill)
            else:
                cc = LineCollection(data, **kwlines)
            ax.add_collection(cc)
            oo.append(cc)
            if isbm:
                m.set_axes_limits(ax=ax)

        # Points
        if points:
            cc = ax.scatter(xx, yy, **kwpoints)
            oo.append(cc)
            if isbm:
                m.set_axes_limits(ax=ax)

        # Special properties
        for key in ['label']:
            if key in kwargs and hasattr(cc, 'set_'+key):
                getattr(cc, 'set_'+key)(kwargs[key])

        # Finalize
        if title:
            ax.set_title(title)
        if autoscale:
            ax.autoscale_view()

        if show: P.show()

        return oo

class XYZ(object):
    """Class to manipulate xyz data (randomly spaced)

    - **xyz**: It can be either

        - a .xyz ascii file, or a netcdf/grd file with variables ``x``, ``y`` and ``z``,
        - a (x,y,z) tuple,
        - a (3, npts) array,
        - another XYZ instance.

    - *long_name*: Long name
    - *units* Units
    - *tranform*: It can be either

        - a factor applied to z at initialisation
        - a fonction that takes z as the only argument to filter its data.

    - *exc*: Polygons to exclude data (see :meth:`exclude`).
             Several polygons must be passed as a tuple (poly1, poly2, ...).
    - *sel*: Polygons to select data (see :meth:`select`).
             Several polygons must be passed as a tuple (poly1, poly2, ...).
    - load_<keywords>: keywords are passed to :func:`numpy.loadtxt`
    - rsamp_<keywords>: keywords are passed to :func:`rsamp`
    - Other keywords are set as atrributes.


    Slicing:

    - len(obj): number of xyz data
    - obj[1:3]: [(x0,y0,z0),(x1,y1,z1)]

    Operations :

    >>> xyz += 2
    >>> xyz3 = xyz1 + xyz2/2. # concatenation
    """

    def __init__(self, xyz, m=None,  units=None, long_name=None, transp=True,
        trans=False, magnet=0, rsamp=0, id=None, **kwargs):
        # Load data
        self._selections = []
        self._exclusions = []
        if isinstance(xyz, XYZ):
            x = xyz._x.copy()
            y = xyz._y.copy()
            z = xyz._z.copy()
            if units is None: units = xyz.units
            if long_name is None: long_name = xyz.long_name
#           if m is None: m = XYZ._m
            self._selections = copy.deepcopy(xyz._selections)
            self._exclusions = copy.deepcopy(xyz._exclusions)
        elif hasattr(xyz, 'xyz'):
            xyz = xyz.xyz
        elif isinstance(xyz, (tuple, N.ndarray)):
            # Direct
            x, y, z = xyz
        elif os.path.exists(xyz):
            # Read from file
            if xyz.endswith('.nc') or xyz.endswith('.grd'):
                # Netcdf or grid file
                f = cdms2.open(xyz)
                x = f('x').filled()
                y = f('y').filled()
                z = f('z').filled()
                f.close()
            else:
                # Ascii file
                data = N.loadtxt(xyz, **kwfilter(kwargs, 'load'))
                x = data[:, 0]
                y = data[:, 1]
                z = data[:, 2]
        else:
            raise TypeError, 'xyz must be either a .xyz file or a tuple of (x, y, z) values'
        # Float64 are needed for polygons
        self._x = N.asarray(x, 'float64')
        self._y = N.asarray(y, 'float64')
        self._z = N.asarray(z)
        if trans is not False:
            if operator.isNumberType(trans):
                self._z *= trans
            elif callable(trans):
                self._z[:] = trans(self._z)
        self._m = None
        self.units = units
        self.long_name = long_name
        self.set_transp(transp)
        self.set_magnet(magnet)
        if rsamp is False: rsamp = 0
        self._rsamp = rsamp
        for att, val in kwargs.items(): setattr(self, att, val)
        self._mask = None
        self._xres_man = self._yres_man = 0
        self._proj = self._proj_(True)
        self.id = None

        # Now update arrays
        self._update_(**kwargs)

    def __str__(self):
        s = 'XYZ:'
        for att in 'long_name', 'units':
            if getattr(self, att) is not None:
                s += '\n %s: %s'%(att, getattr(self, att))
        s += '\n total npts: %i'%self._x.shape[0]
        s += '\n full extension: xmin=%g xmax=%g ymin=%g ymax=%g zmin=%g zmax=%g'%\
            (self._x.min(), self._x.max(), self._y.min(), self._y.max(), self._z.min(), self._z.max())
        s += '\n %i selection polygon(s)'%len(self._selections)
        s += '\n %i exclusion polygon(s)'%len(self._exclusions)
        if len(self._selections) or len(self._exclusions):
            s += '\n filtered extension: xmin=%g xmax=%g ymin=%g ymax=%g zmin=%g zmax=%g'%\
                (self._xmin, self._xmax, self._ymin, self._ymax, self._zmin, self._zmax)
        return s

    def copy(self):
        """Deep copy"""
        self._m = None
        return copy.deepcopy(self)

    def __iadd__(self, other):
        if isinstance(other, (float, int)):
            self._z += other
        else:
            # Consolidate to remove exceptions and selections
            self.consolidate()
            # Load
            other = self._load_(other)
            # Transparency
            if not other.get_transp():
                self.exclude(other)
                self.consolidate()
            # Concatenate arrays
            for l in 'xyz':
                setattr(self, '_%s'%l, N.concatenate((getattr(self, l), getattr(other, l))))
            # Now update everything
            self._update_()
        return self
    def __isub__(self, other):
        assert operator.isNumberType(other), 'Only scalars are allowed for such arithmetic operation'
        self._z -= other
        return self
    def __imul__(self, other):
        assert operator.isNumberType(other), 'Only scalars are allowed for such arithmetic operation'
        self._z *= other
        return self
    def __idiv__(self, other):
        assert operator.isNumberType(other), 'Only scalars are allowed for such arithmetic operation'
        self._z /= other
        return self

    def __add__(self, other):
        newxyz = self.copy()
#       if not isinstance(other, (float, int)):
#           newxyz = self.copy().consolidate()
#           if hasattr(other, '_transp') and not other._transp:
#               # Handle transparency
#               newxyz.exclude(other)
#               newxyz.consolidate()
#           else: # Simply load with no opacity
#               other = self._load_(other)
#       else:
#           newxyz = self.copy()
        newxyz += other
        return newxyz
    def __sub__(self, other):
        newxyz = self._class__(self)
        newxyz -= other
        return newxyz
    def __mul__(self, other):
        newxyz = self._class__(self)
        newxyz *= other
        return newxyz
    def __div__(self, other):
        newxyz = self._class__(self)
        newxyz /= other
        return newxyz

    def _load_(self, d):
        # Load for adding
        # - already a XYZ
        if isinstance(d, self.__class__): return d
        # - XYZ via .xyz() (like shapes/shorelines or mergers)
        if hasattr(d, 'xyz'):
            return d.xyz
        # - load it (from file or array)
        return self.__class__(d)

    def _update_(self,**kwargs):

        # Mask (1==good)
        del self._mask
        npt = len(self._x)
        self._mask = N.ones(npt, '?')
        ii = N.arange(npt, dtype='i')

        # Radius underspampling
        if self._rsamp is False: self._rsamp = 0
        if not hasattr(self, '_rsamp_old'):
            self._rsamp_old = self._rsamp
            self._rsamp_mask = None
        if self._rsamp==0 and self._rsamp_mask is not None:
            del self._rsamp_mask
        elif (self._rsamp!=0 and self._rsamp_mask is None) or self._rsamp != self._rsamp_old:
            del self._rsamp_mask
            kwrsamp=kwfilter(kwargs, 'rsamp')
            self._rsamp_mask = rsamp(self._x, self._y, abs(self._rsamp), proj=self._rsamp<0, getmask=True,**kwrsamp)
        if self._rsamp_mask is not None:
            self._mask &= self._rsamp_mask
        self._rsamp = self._rsamp_old

        # Selections
        if len(self._selections):
            smask = N.zeros(npt, '?')
            rectmask = N.zeros(npt, '?')
            for pl in self._selections:

                # Mask is good only within N/S/W/E limits of the polygon
                rectmask[:] = \
                    (self._x>=pl.boundary[:, 0].min()) & (self._x<=pl.boundary[:, 0].max()) & \
                    (self._y>=pl.boundary[:, 1].min()) & (self._y<=pl.boundary[:, 1].max())

                # Mask is good only within the polygon
                for i in ii[rectmask]:
                    smask[i] |= Point((self._x[i], self._y[i])).within(pl)
            del rectmask
            self._mask &= smask
            del smask

        # Exclusions
        for exc in self._exclusions:

            # Check if inclusions and magnet zone
            if isinstance(exc, tuple):
                exc, incs, magnet = exc
            else:
                incs = []
                magnet = None

            # We check only good points within the N/S/W/E limits of the polygon
            good = self._mask & \
                ((self._x>=exc.boundary[:, 0].min()) & (self._x<=exc.boundary[:, 0].max()) & \
                 (self._y>=exc.boundary[:, 1].min()) & (self._y<=exc.boundary[:, 1].max()))

            # Mask is bad within the polygon
            for i in ii[good]:
                out = not Point((self._x[i], self._y[i])).within(exc)

                # Check inclusions and magnet zone
                if not out:

                    # Inclusions = exclusions of exclusions
                    for inc in incs:
                        out = Point((self._x[i], self._y[i])).within(inc)
                        if out: break

                    # Magnet zone
                    if magnet != None:
                        radius, xmgt, ymgt = magnet
                        dst = N.sqrt((xmgt-self._x[i])**2+(ymgt-self._y[j])**2)
                        out = dst.min() > radius
                        del dst

                # Apply
                self._mask[i] = out
            del good
        del ii

        # Limits
        x = self.x
        y = self.y
        z = self.z
        if len(x):
            self._xmin = x.min()
            self._xmax = x.max()
            self._ymin = y.min()
            self._ymax = y.max()
            self._zmin = z.min()
            self._zmax = z.max()
            self._xmean = x.mean()
            self._ymean = y.mean()
        else:
            warn('No more points after exclusion')
            self._xmin = self._xmax = self._ymin = self._ymax = self._zmin = self._zmax = self._xmean = self._ymean = None

        # Resolution
        if len(x):
            self._xres_mauto, self._yres_mauto = self.resol(deg=False)
            self._xres_dauto, self._yres_dauto = self.resol(deg=True)
        else:
            self._xres_mauto = self._yres_mauto = self._xres_dauto = self._yres_dauto = 0
#       del self._m
#       self._m = None

    def consolidate(self):
        """Apply radius undersampling and all exclusions and selections to data and reset them"""
        x = self.x
        y = self.y
        z = self.z
        self._x = x
        self._y = y
        self._z = z
        self.reset_rsamp()
        self.reset_selections()
        self.reset_exclusions()
        return self

    def mask(self):
        """Get the current mask due to exclusion and selection polygons

        .. seealso::

            :meth:`exclude` :meth:`select`
        """
        return ~self._mask

    def _clip_mask_(self, zone, inverse):
        zone = polygons([zone])[0]
        good = self._mask.copy()
        good = good & \
              (self._x>=zone.boundary[:, 0].min()) & (self._x<=zone.boundary[:, 0].max()) & \
              (self._y>=zone.boundary[:, 1].min()) & (self._y<=zone.boundary[:, 1].max())
        ii = N.arange(len(good), dtype='i')
        for i in ii.compress(good):
            good[i] = Point((self._x[i], self._y[i])).within(zone)
        del ii
        if inverse:
            good = ~good
        return good

    def clip(self, zone=None, margin=None, inverse=False, mask=False, id=None, **kwargs):
        """Geographical selection of part of the data

        - *zone*: (xmin,ymin,xmax,ymax) or a float/int a complex polygon (see :func:`~vacumm.misc.grid.masking.polygons`).
        - *margin*: Margin around ``zone`` relative to the resolution (see :meth:`resol`)
        - *inverse*: Inverse the selection.
        - *mask*: ``zone`` must be interpreted as a mask
        """
        if zone is None or zone == (None, )*4:
            return self
        if margin is not None and isinstance(zone, tuple) and len(zone)==4:
            zone = list(zone)
            if zone[0] is None: zone[0] = self.xmin
            if zone[2] is None: zone[2] = self.xmax
            if zone[1] is None: zone[1] = self.ymin
            if zone[3] is None: zone[3] = self.ymax
            xres, yres = self.get_res()
            xmin = zone[0]-margin*xres
            xmax = zone[2]+margin*xres
            ymin = zone[1]-margin*yres
            ymax = zone[3]+margin*yres
            zone = (xmin, ymin, xmax, ymax)
        if mask is False:
            mask = self._clip_mask_(zone, inverse)
        x = self._x[mask]
        y = self._y[mask]
        z = self._z[mask]
        kwargs.setdefault('units', self.units)
        kwargs.setdefault('long_name', self.long_name)
        kwargs.setdefault('id', id)
        xyz = self.__class__((x, y, z), **kwargs)
#       xyz.include(*self._inclusions)
        xyz.exclude(*self._exclusions)
        return xyz

    def zone(self, poly=False, mask=True):
        """Get xmin,ymin,xmax,ymax

        - *poly*: if True, return zone as a Polygon instance"""
        zone = self.get_xmin(mask), self.get_ymin(mask), self.get_xmax(mask), self.get_ymax(mask)
        if poly: return polygons([zone])[0]
        return zone

#   def get_map(self):
#       """Return the map instance or None"""
#       return self._m


    def set_transp(self, transp):
        """Set :attr:`transp`

        .. note::

            Useful only for mixing :class:`~vacumm.misc.io.XYZ` instances"""
        self._transp = transp
    def get_transp(self):
        """Get :attr:`transp`

        .. note::

            Useful only for mixing :class:`~vacumm.misc.io.XYZ` instances"""
        return self._transp
    transp = property(get_transp, set_transp, doc="Transparency boolean attribute")

    def set_magnet(self, magnet):
        """Set the magnet integer attribute.
        If set to ``0``, no magnet effect.

        .. note::

            Useful only for mixing :class:`~vacumm.misc.io.XYZ` instances"""
        if magnet is False: magnet = 0
        self._magnet = magnet
    def get_magnet(self):
        """Get the magnet integer attribute

        .. note::

            Useful only for mixing :class:`~vacumm.misc.io.XYZ` instances"""
        return self._magnet
    magnet = property(get_magnet, set_magnet, doc="Magnet integer attribute")

    def set_rsamp(self, rsamp):
        """Set the radius sampling :attr:`rsamp`
        If set to ``0``, no sampling."""
        if rsamp is False: rsamp = 0
        self._rsamp = rsamp
        if not hasattr(self, '_rsamp_old') or rsamp != self._rsamp_old:
            self._update_()
    def get_rsamp(self):
        """Get the radius sampling :attr:`rsamp`"""
        return self._rsamp
    def reset_rsamp(self):
        """Reset :attr:`rsamp` without affecting data """
        self._rsamp = self._rsamp_old = 0
        if hasattr(self, '_rsamp_mask'): del self._rsamp_mask
        self._rsamp_mask = None
    del_rsamp = reset_rsamp
    rsmap = property(get_rsamp, set_rsamp, del_rsamp, "Radius of unsersampling")

    def tocfg(self, cfg, section, param=None):
        """Dump one or all parameters as options to a cfg section

        - **cfg**: ConfigParser object
        - **section**: Section of cfg
        - *param*: A single or a list of parameter names
        """
        # List of params
        allowed = ['xmin', 'xmax', 'ymin', 'ymax','zmin','zmax', 'long_name', 'units', 'transp',
                   'xres','yres','exclusions', 'selections']
        if param is None:
            param = allowed
        elif isinstance(param, str):
            param = [param]
        # Get string
        for param in param:
          if param.endswith('res'):
             val = getattr(self, '_'+param+'_mauto')
          elif param.endswith('sions'):
             # selections, exclusions
             val = []
             for var in getattr(self, '_'+param):
                 if isinstance(var, tuple):
                    # exclusions and inclusions
                    exc, incs = var
                    if len(incs) == 0:
                       # no inclusions
                       val.append(exc.get_coords().tolist())
                    else:
                       polys = []
                       for inc in incs:
                          polys.append(inc.get_coords().tolist())
                          val.append((exc, polys))
                 else:
                     # normal case
                     val.append(var.get_coords().tolist())
          elif param in ['units', 'long_name']:
              val = getattr(self, param)
          else:
              val = getattr(self, '_'+param)
          # Check section
          if not cfg.has_section(section):
              cfg.add_section(section)
          # Dump
          cfg.set(section, param, str(val))


    def get_xmin(self, mask=True):
        if mask is True or mask == 'masked': return self._xmin
        return self.get_x(mask).min()
    xmin = property(get_xmin, doc="X min")
    def get_xmax(self, mask=True):
        if mask is True or mask == 'masked': return self._xmax
        return self.get_x(mask).max()
    xmax = property(get_xmax, doc="X max")
    def get_ymin(self, mask=True):
        if mask is True or mask == 'masked': return self._ymin
        return self.get_y(mask).min()
    ymin = property(get_ymin, doc="Y min")
    def get_ymax(self, mask=True):
        if mask is True or mask == 'masked': return self._ymax
        return self.get_y(mask).max()
    ymax = property(get_ymax, doc="Y max")
    def get_zmin(self, mask=True):
        if mask is True or mask == 'masked': return self._zmin
        return self.get_z(mask).min()
    zmin = property(get_zmin, doc="Z min")
    def get_zmax(self, mask=True):
        if mask is True or mask == 'masked': return self._zmax
        return self.get_z(mask).max()
    zmax = property(get_zmax, doc="Z max")

    def exclude(self, *zones):
        """Add one or more zones where data are not used.

        A zone can be :

        - an argument to :func:`~vacumm.misc.grid.masking.polygons` to get a :class:`_geoslib.Polygon` instance,
        - another :class:XYZ` instance from which the convex hull (see :meth:`hull`) is used as a delimiting area

        :Usage:

        >>> xyz.exclude([[-8,43],[-5.5,43],[-6,45.]],[[-10,45],[-7,47],[-10,49.]])
        >>> xyz.exclude(polygon1,polygon2)
        >>> xyz.exclude(xyz1,[-5,42,-3,48.])

        .. seealso::

            :meth:`select` :meth:`exclusions`
        """
        zones = list(zones)
        for i, zone in enumerate(zones):
            if isinstance(zone, XYZ):
                zone = zone.shadows()
            if isinstance(zone, tuple): # inclusions
                exc = polygons([zone[0]])[0]
                if len(zone) == 1 or len(zone[1])==0:
                    zone = exc
                else:
                    zone = exc, polygons(zone[1]), zone[2]
            else:
                zone = polygons([zone])[0]
            zones[i] = zone
        self._exclusions.extend(zones)
        self._update_()
    def select(self, *zones):
        """Add one or more zone (polygons) where only these data are used

        A zone is an argument to :func:`~vacumm.misc.grid.masking.polygons` to get a :class:`_geoslib.Polygon` instance.

        :Usage:

        >>> xyz.select([[-8,43],[-5.5,43],[-6,45.]],[[-10,45],[-7,47],[-10,49.]])
        >>> xyz.select(polygon1,polygon2)

        .. seealso::

            :meth:`exclude`  :meth:`selections`
        """
        self._selections.extend(polygons(list(zones)))
        self._update_()
    def reset_selections(self):
        """Remove all selections"""
        del self._selections
        self._selections = []
        self._update_()
    def reset_exclusions(self):
        """Remove all exclusions"""
        del self._exclusions
        self._exclusions = []
        self._update_()
    def selections(self):
        """Get all selection polygons as a tuple"""
        return self._selections
    def exclusions(self):
        """Get all exclusion polygons as a tuple"""
        return self._exclusions


    def _filter_(self, var, mask):
        """
        - var can be a list
        - mask can be
            - 'valid', True
            - False, None
            - 'revert', 'masked'
            - mask array
        """
        if mask == 'valid': mask = True
        assert self._mask is not None
        if mask is False :
            return var
        if isinstance(mask, str) and (mask.startswith('rever') or mask.startswith('inver') or
            mask.startswith('mask')):
            mask = ~self._mask
        else:
            mask = self._mask
        if isinstance(var, list):
            return [v.compress(mask) for v in var]
#       if isinstance(var, list):
#           return [var[i] for i, m in enumerate(self._mask) if m]
        return var.compress(mask)
    def get_x(self, mask=True):
        """Get valid X positions"""
        return self._filter_(self._x, mask)
    x = property(get_x, doc='Valid X positions')
    def get_y(self, mask=True):
        """Get valid Y positions"""
        return self._filter_(self._y, mask)
    y = property(get_y, doc='Valid Y positions')
    def get_z(self, mask=True):
        """Get valid Z values"""
        return self._filter_(self._z, mask)
    z = property(get_z, doc='Valid Z values')
    def get_xy(self, mask=True):
        """Return coordinates as a (2, npts) array :attr:`xy`

        - ``xy()[0]``: X
        - ``xy()[1]``: Y
        """
        return N.asarray([self.get_x(mask), self.get_y(mask)])
    xy = property(get_xy, doc='Coordinates as a (2, npts) array')
    def get_xyz(self, mask=True, split=False):
        """Return coordinates and data as a (3, npts) array :attr:`xyz`

        - ``xy()[0]``: X
        - ``xy()[1]``: Y
        - ``xy()[2]``: Z
        """
        if split: return self.get_x(mask), self.get_y(mask), self.get_z(mask)
        return N.asarray([self.get_x(mask), self.get_y(mask), self.get_z(mask)])
    xyz = property(get_xyz, doc='Coordinates and data as a (3, npts) array')

    def astuple(self, mask=True):
        """Shortcut to ``xyz(split=True)`` (see :meth:`xyz`)"""
        return self.get_xyz(mask=mask, split=True)


    def __len__(self):
        return self._mask.sum()
    def __getitem__(self, key):
        x = self.x[key]
        y = self.y[key]
        z = self.z[key]
        if isinstance(key, int):
            return x, y, z
        return zip(x, y, z)

    def interp(self, xyo, xyz=False, **kwargs):
        """Interpolate to (xo,yo) positions using :class:`nat.Natgrid`

        :Params:

            - **xo**: Output X
            - **yo**: Output Y
            - *xyz*: If True, return a :class:`XYZ` instance instead of a :mod:`numpy` array
            - **interp_<param>**, optional: ``<param>`` is passed to the
              :func:`~vacumm.misc.grid.regridding.xy2xy` interpolation routine.
            - Other params are passed to XYZ initialization for the output dataset.

        :Returns: An XYZ instance
        """
        # FIXME: Natgrid still required by this module ??
        # from nat import Natgrid
        if isinstance(xyo, (tuple, N.ndarray)):
            xo, yo = xyo
        elif hasattr(xyo, 'x'):
            xo = xyo.x
            yo = xyo.y
        elif hasattr(xyo, 'xy'):
            xo, yo = xyo.xy
        elif hasattr(xyo, 'xyz'):
            xo, yo, tmp = xyo.xyz
        else:
            raise TypeError, 'Wrong input type'
        kwinterp = kwfilter(kwargs, 'interp_')
        from grid.regridding import xy2xy
        zo = xy2xy(self.x, self.y, self.z, xo, yo, **kwinterp)
        if not xyz: return zo
        kwargs.setdefault('units', self.units)
        kwargs.setdefault('long_name', self.long_name)
        return self.__class__((xo, yo, zo), **kwargs)

    def hull(self, out='xy', mask=True):
        """Return the convex hull

        :Returns: Depends on ``out``

        - ``"xy"``: (xhull, yhull)
        - ``"ind"``: indices of points
        - ``"poly"``: :class:`_geoslib.Polygon` instance
        """
        return convex_hull((self.get_x(mask), self.get_y(mask)), poly=out=='poly')

    def shadows(self):
        """Get the polygons defining the 'shadow' of this dataset.

        It consists of a tuple of two elements:

            - the convex hull as a polygon,
            - a list of exclusion polygons that intersect the convex hull.

        Therefore, a point in the shadow must be inside the convex hull polygon,
        and outside the exclusion polygons.

        :Returns: (hull_poly, [exclusion_poly1,...])
        """
        hull_poly = self.hull(out='poly')
        exc_polys = []
        for exc in self.exclusions():
            if isinstance(exc, tuple): exc = exc[0]
            if hull_poly.intersects(exc):
                exc_polys.append(exc)
            #FIXME: inclusions
        magnet = None
        if self.get_magnet():
            if self.get_magnet() < 0:
                radius = -self.get_magnet()
            else:
                xres, yres = self.get_res(deg=True)
                radius = N.sqrt((xres**2+yres**2)/2.)*self.get_magnet()
            magnet = radius, self.x, self.y
        return hull_poly, exc_polys, magnet


    def contains(self, x, y):
        """Check if one or several points are within a the convex hull

        - **x,y**: X,Y positions as floats or lists or an :mod:`numpy` arrays.

        .. seealso:

            :meth:`hull`
        """
        hull = self.hull(out='poly')
        try:
            len(x)
            x = N.asarray(x)
            y = N.asarray(y)
            good = (x>=self.xmin) & (x<=self.xmax) & \
                  (y>=self.ymin) & (y<=self.ymax)
            xg = x.compress(good)
            yg = y.compress(good)
            out = []
            for xp, yp in zip(xg, yg):
                out.append(Point(xp, yp).within(hull))
            if isinstance(x, list): return out
            return N.array(out)
        except:
            if (x<self.xmin)|(x>self.xmax)|(y<self.ymin)|(y>self.ymax):
                return False
            return Point(x, y).within(hull)


    def resol(self, convex_hull_method='delaunay', exc=[], deg=False):
        """Return the mean resolution.

        Algorithm: Median distances between facets of triangles

        :Returns: (xres,yres)

        """
        # Coordinates
        x = self.x
        y = self.y
        if not deg:
            x, y = self._proj_(True)(x, y)

#       if method.startswith('tri'): # Using the median length of triangles
        from scipy.spatial import Delaunay

        coord=N.vstack((x,y)).T
        t=Delaunay(coord)
        distances = []
        for p1,p2,p3 in t.vertices:
           distances.append(N.sqrt((x[p1]-x[p2])**2+(y[p1]-y[p2])**2))
           distances.append(N.sqrt((x[p1]-x[p3])**2+(y[p1]-y[p3])**2))
           distances.append(N.sqrt((x[p2]-x[p3])**2+(y[p2]-y[p3])**2))

        xres = yres = N.median(distances)
        del t, distances

#       else: # Using the convex hull
#
#           # Convex hull polygon
#           hull = convex_hull((x, y), poly=True, method=convex_hull_method)
#
#           # Area
#           # - max = convex hull
#           area = hull.area()
#           # - substract intersections with exclusions
#           #FIXME: xyz.resol: non overlapping exclusions+inclusions
#           for e in exc+self.exclusions():
#               if isinstance(e, tuple):
#                   e, incs = e
#               else:
#                   incs = []
#               if hull.intersects(e):
#                   for i in hull.intersection(e):
#                       area -= i.area()
#
#           # Area and mean resolution
#           xres = yres = N.sqrt(area/len(x))

        return xres, yres

    def set_res(self, xres, yres=None):
        """Set the resolution of the dataset

        If ``yres`` is not, it is set to ``xres``.
        When a value is **negative**, it is supposed to be in **meters** (not in degrees)"""
        if yres is None: yres = xres
        self._xres_man = xres
        self._yres_man = yres

    def get_res(self, deg=False, auto=None):
        """Get the mean X and Y resolutions in meters or degrees"""
        # Get manual or auto resolutions
        if self._xres_man and auto is not True:
            xres = self._xres_man
        elif deg:
            xres = self._xres_dauto
        else:
            xres = self._xres_mauto
        if self._yres_man and auto is not True:
            yres = self._yres_man
        elif deg:
            yres = self._yres_dauto
        else:
            yres = self._yres_mauto
        # Conversions
        return self._get_res_(xres, yres, deg)

    def _get_res_(self, xres, yres, degres):
        if degres: # need degrees
            if xres<0: xres = -m2deg(xres, self._ymean)
            if yres<0: yres = -m2deg(yres)
            return xres, yres
        else: # need meters
            if xres>0: xres = deg2m(xres, self._ymean)
            if yres>0: yres = deg2m(yres)
            return xres, yres


    def _proj_(self, proj):
        """Get a proper projection or None"""
        if proj is None: proj = False
        if not proj: return
        if not callable(proj):
            from vacumm.misc.grid.basemap import get_proj
            proj = get_proj((self._x, self._y))
        return proj

    def get_grid(self, res=None, xmin=None, xmax=None, ymin=None, ymax=None, relres=.5, degres=False, id='xyz_grid'):
        """Generate a rectangular grid based on x/y positions and resolution

        - *res*: Resolution. It can be:

                - a float where then ``xres=yres=res``
                - a tuple as ``(xres,yres)``
                - else it is guessed using :meth:`get_res` (and maybe :meth:`resol`)` and multiplied by ``relres``

        - *relres*: Relative resolution factor applied to ``res`` when resolution is guessed (``res=None``)
        - *degres*: When ``res`` is explicitly given, it interpreted as degrees is ``degres`` is True.
        - *xmin,xmax,ymin,ymax*: Bounds of the grid. If not specified, bounds of the dataset are used (see :meth:`xmin`, etc).

        .. note::

            Resolutions are adjusted when they are not mutiple of grid extensions (slightly decreased).
            Therefore, extensions of the grid are always preserved.

        .. seealso::

            :meth:`resol`, :meth:`togrid`
        """
        # Resolution
#       proj = self._proj_(proj)
        if res is None: # Auto
            # Meters
#           dx, dy = tuple([r*relres for r in self.get_res(deg=False, auto=True)])
            dx, dy = tuple([r*relres for r in self.get_res(deg=True, auto=True)])
#           # To degrees
#           dx, dy = self._get_res_(-dx, -dy, True)
        else: # Manual
            if isinstance(res, tuple):
                xres, yres = res
            else: # Same for X and Y
                xres = yres = res
            if not degres: # Given in meters
                xres = -xres
                yres = -yres
            # Now in degrees
            dx, dy = self._get_res_(xres, yres, True)

        # Bounds
        if xmin is None: xmin = self.xmin
        if xmax is None: xmax = self.xmax
        if ymin is None: ymin = self.ymin
        if ymax is None: ymax = self.ymax

        # Rectifications
        xratio = (xmax-xmin)/dx
        yratio = (ymax-ymin)/dy
        dx += dx*(xratio%1)/int(xratio)
        dy += dy*(yratio%1)/int(yratio)

        # Grid
        grid = create_grid(N.arange(xmin, xmax+dx/2., dx), N.arange(ymin, ymax+dy/2., dy))
        grid.id = id
        return grid
    grid = property(get_grid, doc="Rectangular grid based on x/y positions and resolution")

    def togrid(self, grid=None, mask=False, cgrid=False,  **kwargs):
        """Interpolate to a regular grid

        - **grid**: The output grid. It can be either:

            - a (x,y) tuple or a grid or a :mod:`MV2` variable with a grid,
            - ``None``, thus guessed using :meth:`grid`

        - *mask*: It can be either:

            - ``None``, ``False`` or ``MV2.nomask``: no masking
            - an array: this mask array is directly applied
            - a :class:`Shapes` instance (or :class:`~vacumm.bathy.shorelines.ShoreLine`) or a single char GSHHS resolution (and optionally 's' for Histolitt)
            - a callable fonction so that ``mask = thisfunc(mask, **kwmask)``
            - a float: data with this value are masked

        - *mask_<param>*: <param> is passed to :func:`~vacumm.misc.grid.masking.polygon_mask` for evaluation of mask thanks to the polygons.
        - *grid_<param>*: <param> is passed to :func:`grid`.
        - *cgrid*: If ``True``, returns bathy at U- and V-points, else at T-points
        - Other keyparam are passed to :func:`~vacumm.misc.grid.regridding.griddata` for regridding.

        Return: ``(Zx,Zy)`` OR ``Z`` depending on cgrid.

        .. seealso:

            :func:`~vacumm.misc.grid.masking.polygon_mask` :class:`~vacumm.misc.grid.basemap.GSHHS_BM`
            :class:`Shapes`  :class:`~vacumm.bathy.shorelines.ShoreLine`
        """
        # Grid
        kwgrid = kwfilter(kwargs, 'grid_')
        if grid is None:
            grid = self.get_grid(**kwgrid)

        # Interpolation
        kwmask = kwfilter(kwargs, 'mask_')
        kwargs.setdefault('method', 'linear')
        kwargs.setdefault('ext', True)
        var = griddata(self.x, self.y, self.z, grid, mask=None, cgrid=cgrid, **kwargs)
        if not cgrid: var = var,

        # Mask from polygons
        if isinstance(mask, (Shapes, str)):

            # Clip for polygons
            clip = kwmask.pop('clip', .1)
            if isinstance(clip, float):
                xx, yy = get_xy(grid)
                xmin = xx[:].min() ; xmax = xx[:].max()
                ymin = yy[:].min() ; ymax = yy[:].max()
                dx = (xmax-xmin)*clip
                dy = (ymax-ymin)*clip
                clip = [xmin-dx, ymin-dy, xmax+dx, ymax+dy]

            # Get polygons
            if isinstance(mask, Shapes): # Direct
                shapes = mask.clip(clip)
            else: # GSHHS
                from grid.basemap import GSHHS_BM
                shapes = GSHHS_BM(mask, clip=clip)

            # Create mask
            mask = polygon_mask(grid, shapes.get_shapes(), **kwmask)

        # Masking
        if mask is not None and mask is not False and mask is not MV2.nomask:
            for vv in var:
                if hasattr(mask, 'ndim'): # array
                    vv[:] = MV.masked_where(mask, vv, copy=0)
                elif callable(mask): # function
                    vv[:] = mask(vv, **kwmask)
                else: # value
                    vv[:] = MV.masked_values(vv, mask, copy=0)

        # Attributes
        for vv in var:
            if self.units is not None:
                vv.units = self.units
            if self.long_name is not None:
                vv.long_name = self.long_name
        if cgrid: return var
        return var[0]

    def toxy(self, xo, yo, mask=None, outtype='tuple'):
        """Interpolate on random points using :func:`~vacumm.misc.grid.regridding.xy2xy`

        - **xo,yo**: Output positions
        - *mask*: It can be either:

            - ``None``, ``False`` or ``MV2.nomask``: no masking
            - a :class:`Shapes` instance (or :class:`~vacumm.bathy.shorelines.ShoreLine`) or a single char GSHHS resolution (and optionally 's' for Histolitt)
        - *outtype*: Define output type

            - ``"tuple"``: as a tuple (x, y, z)
            - ``"xyz"``: as xyz block
            - ``"XYZ"``: as an :class:`XYZ` (or subclass) instance
        """

        # Interpolate
        zo = xy2xy(self.x, self.y, self.z, xo, yo, **kwargs)


        # Mask from polygons
        if isinstance(mask, (Shapes, str)):

            # Clip for polygons
            clip = kwmask.pop('clip', .1)
            if isinstance(clip, float):
                dx = (xo.max()-xo.min())*clip
                dy = (yo.max()-yo.min())*clip
                clip = [xo.min()-dx, yo.min()-dy, xo.max()+dx, yo.max()+dy]

            # Get polygons
            if isinstance(mask, Shapes): # Direct
                shapes = mask.clip(clip)
            else: # GSHHS
                from grid.basemap import GSHHS_BM
                shapes = GSHHS_BM(mask, clip=clip)

            # Maskit
            xo, yo, zo = polygon_select(xo, yo, shapes.get_shapes(), zz=zo)

        # Out
        if outype=='xyz':
            return N.asarray([xo, yo, zo])
        if outtype=='XYZ':
            return self.__class__((xo, yo, zo))
        if callable(outtype):
            return outtype([xo, yo, zo])
        return xo, yo, zo



    def plot(self, size=5., color=None, alpha=1., masked_alpha=.3,
        masked_size=None, linewidth=0., show=True, savefig=None,
        savefigs=None, m=None, colorbar=True, title=None,
        units=None,  cmap=None, mode='valid', zorder=100,
        masked_zorder=50, margin=2,
        xmin=None, xmax=None, ymin=None, ymax=None,
        xres=None, yres=None, **kwargs):
        """Scatter plot of bathymetry points

        :Params:

            - **mode**, optional: 'valid', 'masked' or 'both'.
            - **size**, optional: Size of markers.
            - **color**, optional: Color of markers.
            - **alpha**, optional: Alpha transparency of markers.
            - **zorder**, optional: zorder of markers.
            - **m**, optional: Use this :class:`~mpl_toolkits.basemap.Basemap` instance
              to plot the points.
            - **masked_size**, optional: Size of masked markers.
            - **masked_alpha**, optional: Alpha transparency of masked markers.
            - **masked_zorder**, optional: zorder of masked markers.
        """

        # Params
        kwm = kwfilter(kwargs, 'map')
        kwhull = kwfilter(kwargs, 'hull')
        kwcb = kwfilter(kwargs, 'colorbar')
        kwplot = dict(linewidth=linewidth, cmap=get_cmap(cmap))
        kwplot.update(kwargs)
        pts = None

        # Limits
        if m is True:
            if xmin is None: xmin = self.xmin
            if xmax is None: xmax = self.xmax
            if ymin is None: ymin = self.ymin
            if ymax is None: ymax = self.ymax
        if margin != 0:
            if xres is None or yres is None:
                xxres, yyres = self.resol(deg=True)
            if xres is None: xres = xxres
            if yres is None: yres = yyres
            if xmin is not None: xmin -= margin*xres
            if xmax is not None: xmax += margin*xres
            if ymin is not None: ymin -= margin*yres
            if ymax is not None: ymax += margin*yres

        # Map
        if m is True:
            kwm.setdefault('lon', (xmin, xmax))
            kwm.setdefault('lat', (ymin, ymax))
            kwm['projection'] = 'merc'
            m = map2(show=False, savefig=None, **kwm)
        if m:
            G = m.map if hasattr(m, 'map') else m
            self._m = G
        else:
            G = P

        # Max values
        if mode == 'both':
            zmin = self._z.min()
            if not kwplot.has_key('vmin'): kwplot['vmin'] = zmin
            zmax = self._z.max()
            if not kwplot.has_key('vmax'): kwplot['vmax'] = zmax

        # Valid points
        if mode != 'masked' or mode == 'both':
            # Points
            pts = self._plot_('valid', G, m, size, color, alpha=alpha, zorder=zorder,
                label=self.long_name, **kwplot)
#           # Convex hull
#           if hull:
#               xhull, yhull = self.hull(out='xy')
#               if callable(m):
#                   xhull, yhull = m(xhull, yhull)
#               xhull = xhull.tolist()+[xhull[0]]
#               yhull = yhull.tolist()+[yhull[0]]
#               kwhull.setdefault('alpha', masked_alpha)
#               kwhull.setdefault('zorder', zorder)
#               if color is None:
#                   hullcolor = '#888888'
#               else:
#                   hullcolor = color
#               kwhull.setdefault('color', hullcolor)
#               G.plot(xhull, yhull, '-', **kwhull)
        # Masked points
        if mode == 'masked' or mode == 'both':
            if mode == 'masked':
                masked_alpha = alpha
                masked_size = size
                masked_zorder = zorder
            elif masked_size is None: masked_size = size
            p = self._plot_('masked', G, m, masked_size, color, alpha=masked_alpha, zorder=masked_zorder, **kwplot)
            if pts is None: pts = p

        # Limits
        if m:
            G.set_axes_limits(P.gca())
        else:
            if xmin is not None: P.xlim(xmin=xmin)
            if xmax is not None: P.xlim(xmax=xmax)
            if ymin is not None: P.ylim(ymin=ymin)
            if ymax is not None: P.ylim(ymax=ymax)

        # Decorations
        if title is None and self.long_name is not None:
            title = self.long_name
        if title is not None:
            if self.long_name is None:
                self.long_name = title
            P.title(title)
#           for pt in pts:
#               if pt is not None:
#                   pt.set_label(title)
#                   break
        if units is None and self.units is not None:
            units = self.units
        if units is not None:
            if self.units is None:
                self.units = units
        if colorbar:
            _colorbar_(pts, units=units, **kwcb)
        if savefig:
            P.savefig(savefig)
        elif savefigs:
            Savefigs(savefigs)
        if show:
            P.show()
        return pts

    def _plot_(self, mask, G, proj, size, color, **kwargs):
        """Generic plot"""
        x, y = self.get_x(mask), self.get_y(mask)
        if not len(x): return None
        if callable(proj): x, y,  = proj(x, y)
        if color is None: color = self.get_z(mask)
        kwargs.setdefault('linewidth', 0)
        if kwargs.get('marker', None) is None:
            kwargs['marker'] = 'o'
        return  G.scatter(x, y, s=size, c=color, **kwargs)

    def save(self, xyzfile, **kwargs):
        """Save to a file

        - **xyzfile**: Output file name

            - write a netcdf file if it ends with ".nc" or ".grd"
            - write a sinux file if it ends with ".snx"
            - else write an ascii file with 3 columns

        - Other keywords are passed to :func:`numpy.savetxt` for ascii saving
        """
        if xyzfile.endswith('.nc') or xyzfile.endswith('.grd'):
            x = cdms2.createVariable(self.x, id='x')
            y = cdms2.createVariable(self.y, id='y')
            z = cdms2.createVariable(self.z, id='z')
            i = cdms2.createAxis(range(len(x)), id='i')
            f = cdms2.open(xyzfile, 'w')
            for var in x, y, z:
                var.setAxis(0, i)
                f.write(var)
            f.close()
        elif xyzfile.endswith('.snx'):
            write_snx(self.xyz.T, type='point', **kwargs)
        else:
            N.savetxt(xyzfile, self.xyz.T, **kwargs)


class XYZMerger(object):
    """Mix different bathymetries"""
    def __init__(self, *datasets, **kwargs):
        self._xmin = kwargs.pop('xmin', None)
        self._xmax = kwargs.pop('xmax', None)
        self._ymin = kwargs.pop('ymin', None)
        self._ymax = kwargs.pop('ymax', None)
        self._datasets = list(datasets)
        self._XYZ = XYZ
        self.long_name = kwargs.get('long_name', None)
        self.units = kwargs.get('units', None)

    def __str__(self):
        s = 'Merger:'
        for att in 'long_name', 'units':
            if getattr(self, att) is not None:
                s += '\n %s: %s'%(att, getattr(self, att))
        n = len(self)
        if not n:
            return s+'\n no datasets in the merger'
        xyz = self.xyz
        s += '\n extension: xmin=%g xmax=%g ymin=%g ymax=%g zmin=%g zmax=%g'%\
            (xyz.xmin, xyz.xmax, xyz.ymin, xyz.ymax, xyz.zmin, xyz.zmax)
        if n==1:
            s += '\n there is 1 dataset in the merger:'
        else:
            s += '\n there are %i datasets in the merger:'%n
        sep = '\n'+'-'*3
        s += sep
        s += sep.join(['\n'+str(d) for d in self])
        s += sep
        return s


#   def __copy__(self):
#       return self.__class__(xmin=self._xmin, xmax=self._xmax, ymin=self._ymin, ymax=self._ymax, *self._datasets)
    def copy(self):
        return copy.deepcopy()
#       return self.__copy__()

    def _load_(self, d):
        # XYZMerger instance
        if isinstance(d, self.__class__):
        # xyzmerger._load_: pas clair
#           if d is self: return
#           return d.copy()
            return [self._XYZ(dd) for dd in d._datasets]
        # XYZ instance
        if isinstance(d, self._XYZ):
            return d
        # Get XYZ from b
        if hasattr(d, 'xyz'):
            return d.xyz
        # New XYZ instance
        return self._XYZ(d)

    def tolist(self):
        """Return the merger as a list of datasets"""
        return self._datasets

    def ids(self):
        return [d.id for d in self._datasets]

    def append(self, d):
        """Append a dataset to the merger"""
        self += d
    def remove(self, d):
        """Remove a dataset from the merger"""
        self -= d
    def clean(self):
        """Remove all current dataset"""
        del self._datasets
        self._datasets = []

    def __iadd__(self, datasets):
        if not isinstance(datasets, list): datasets = [datasets]
        for d in datasets:
            d = self._load_(d)
            if d is not None and d not in self._datasets:
                self._datasets.append(d)
        return self
    def __isub__(self, d):
        if d not in self._datasets:
            warn('Dataset not in merger')
        self._datasets.remove(d)
        return self

    def __getitem__(self, key):
        return self._datasets[self._key_(key)]
    def __delitem__(self, key):
        del self._datasets[self._key_(key)]
    def __setitem__(self, key, b):
        self._datasets[self._key_(key)] = self._load_(b)
    def __len__(self):
        return len(self._datasets)
    def _key_(self, key):
        if isinstance(key, str):
            key = self.ids().find(key)
        return key

    def get_xyz(self, mask=True, **kwargs):
        """Merge current dataset"""
        assert len(self), 'You must add at least one dataset to the merger'
        xyz = self._XYZ(self._datasets[0], **kwargs)
        if mask is False:
            xyz.reset_exclusions()
            xyz.reset_selections()
        if mask:
            pass
        for d in self[1:]:
            if mask is False:
                d = copy.deepcopy(d)
                d.reset_exclusions()
                d.reset_selections()
            xyz += d
            if mask:
                pass
        xyz.long_name = self.long_name
        xyz.units = self.units
        return xyz
    xyz = property(get_xyz, doc='Coordinates and data as a (3, npts) array')

    def merge(self, **kwargs):
        """Shortcut to :meth:`xyz`"""
        return self.get_xyz(**kwargs)

    def togrid(self, *args, **kwargs):
        """Interpolate merged bathymetries to a grid"""
        return self.xyz.togrid(*args, **kwargs)

    def plot(self, color=None, marker=None, mode='cluster', title='XYZ merger',
        show=True, colorbar=True, savefig=None, savefigs=None, legend=True,
        xmin=None, xmax=None, ymin=None, ymax=None, margin=5, xres=None, yres=None, **kwargs):
        """

        - *alpha*: Alpha transparency:

            - applied to **all** points if ``mode="cluster"``
            - applied to **hidden** points if ``mode="data"``

        - *mode*: Display mode:

            - ``"cluster"``: Points from different datasets have different colors and markers,
                and hidden points are transparent.
            - ``"data"``: Points have the same marker, colors depends on Z value and hidden
                points are masked.

        - *marker*: Define a single or several markers to be used.
        - *legend*: Show a legend if ``mode="cluster"``.
        - *title*: Title of the plot.
        - *m*: :class:`~mpl_toolkits.basemap.Basemap` instance.
        - *m_margin*: Margin for ``m``, relative to the mean resolution (see :meth:`XYZ.resol`)
        - *m_<keywords>*: Keywords are passed to :func:`~vacumm.misc.plot.map`.
        - Extra keywords are passed to :meth:`XYZ.plot`.
        """


        # Limits
        if None in [xmin, xmax, ymin, ymax] or margin != 0.:
            xyz = self.get_xyz(mask=False)
        if xmin is None: xmin = xyz.get_xmin(False)
        if xmax is None: xmax = xyz.get_xmax(False)
        if ymin is None: ymin = xyz.get_ymin(False)
        if ymax is None: ymax = xyz.get_ymax(False)
        if margin != 0. and None in [xres, yres]:
            xxres, yyres = xyz.resol(deg=True)
            if xres is None: xres = xxres
            if yres is None: yres = yyres

        # Arguments to plots
        kwleg = kwfilter(kwargs, 'legend')
        kwsf = kwfilter(kwargs, 'savefig')
        kwsfs = kwfilter(kwargs, 'savefigs')
        kwplot = dict(show=False, colorbar=False, title=False,
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            margin=margin, xres=xres, yres=yres)
        kwplot.update(kwargs)

        # Select mode
        pp = []
        if mode == 'cluster':

            # Colors and symbols
            if color is None:
                color = simple_colors
            color = broadcast(color, len(self))
            if marker is None:
                marker = Markers
            marker = broadcast(marker, len(self))

            # Copy all datasets
            datasets = [copy.deepcopy(d) for d in self._datasets]

            # Loop on datasets
            m = kwplot.pop('m', None)
            for j, i in enumerate(xrange(len(datasets)-1, -1, -1)):

                # Plot current dataset
                pp.append(datasets[i].plot(mode='both', color=color[i], marker=marker[i],
                    zorder=200+i, masked_zorder=50+i, m=m, **kwplot))
                m = datasets[i]._m

                # Mask next datasets with current one if not transparent
                if i != 0 and not datasets[i].get_transp():
                    for k in xrange(i-1, -1, -1):
                        datasets[k].exclude(datasets[i])

            # Legend
            if legend:
                hh = ()
                ll = ()
                for d, p in zip(datasets[::-1], pp):
                    if d.long_name is not None:
                        hh += p,
                        ll += d.long_name,
                legend_alpha = kwleg.pop('alpha', .3)
                kwleg.setdefault('loc', 'best')
                kwleg.setdefault('shadow', False)
                leg = P.legend(hh, ll, **kwleg)
                leg.legendPatch.set_alpha(legend_alpha)
                leg.set_zorder(1000)

        else:

            # Physical case
            pp.append(self.xyz.plot(color=color, marker=marker, **kwplot))
            if colorbar:
                _colorbar_(pp[0])

        # Plot attributes
        if title is None and self.long_name is not None:
            title = self.long_name
        if title is not False:
            P.title(title)
        if savefig:
            P.savefig(savefig, **kwsf)
        if savefigs:
            Savefigs(savefigs, **kwsfs)
        if show:
            P.show()
        return pp



def write_snx(objects, snxfile, type='auto', mode='w', z=99, xfmt='%g', yfmt='%g', zfmt='%g', close=True):
    """Write points, lines or polygons in a sinusX file"""
    # Check depth
    if isinstance(objects, (LineString, Polygon)):
        objects = [objects]
    elif isinstance(objects[0], Point) or \
        (not isinstance(objects[0], (LineString, Polygon)) and not hasattr(objects[0][0], '__len__')):
        objects = [objects]
        if type =='auto' or isinstance(objects[0][0], Point): type = 'point'

    # File
    splitfile = False
    if isinstance(snxfile, file):
        f = snxfile
    elif '%i' not in snxfile:
        f = open(snxfile, mode)
    else:
        splitfile = True
        snxfile = snxfile.replace('%i', '%%0%ii'%int(N.log10(len(objects))))

    # Loop on objects
    for i, oo in enumerate(objects):

        # Extract object
        if isinstance(oo, LineString):
            oo = oo.get_coords()
            type = 'linestring'
        elif isinstance(oo, Polygon):
            oo = oo.get_coords()
            type = 'polygon'
        elif isinstance(oo[0], tuple):
            type = 'point'
        elif isinstance(oo[0], Point):
            ooo = []
            for o in oo:
                ooo.extend(o.get_coords())
            oo = ooo
            type = 'point'

        # Guess type
        if type == 'auto':
            if isinstance(oo[0], float):
                type = 'point'
            else:
                n = N.asarray(oo)
                d = N.sqrt(N.diff(n[:, 0])**2+N.diff(n[:, 1])**2)
                if d.mean() < 3*N.sqrt((n[0, 0]-n[-1, 0])**2+(n[0, 1]-n[-1, 1])**2):
                    type='polygon'
                else:
                    type = 'linestring'
                del d, n

        # Write
        # - splited file
        if isinstance(snxfile, (str, unicode)) and '%' in snxfile:
            f = open(snxfile%i, mode)
        # - header
        if type.startswith('point'): # Points
            f.write("B S\nCN Semis\nCP 0 0\nCP 0\n")
        elif type.startswith('line'): # LineString
            f.write("B N\nCN Niveau\nCP 0 1\nCP 99\nCP 0\n")
        elif type.startswith('poly'): # Polygon
            f.write("B N\nCN Niveau\nCP 1 1\nCP 99\nCP 0\n")
        for o in oo:
            if len(o) == 2:
                zz = z
            else:
                zz = o[2]
            f.write(('%s %s %s'%(xfmt, yfmt, zfmt)+' A\n')%(o[0], o[1], zz))
        if isinstance(snxfile, (str, unicode)) and '%' in snxfile:
            f.close()
    if not f.closed and close:
        f.close()


class ColoredFormatter(logging.Formatter):
    """Log formatter with colors"""
    def __init__(self, msg, full_line=False):
        logging.Formatter.__init__(self, msg)
        self.full_line = full_line
        self.colorize = TermColors().format

    def format(self, record):
        if self.full_line:
            return self.colorize(logging.Formatter.format(self, record), record.levelname)
        record.levelname = self.colorize(record.levelname, record.levelname)
        return logging.Formatter.format(self, record)

class _Redirector_(object):
    def __init__(self, func, prefix=''):
        self.func = func
        self.prefix = prefix
    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.func(self.prefix+line.rstrip())
    def flush(self):
        pass




class Logger(object):
    """Class for logging facilities when subclassing.
    Logging may be sent to the console and/or a log file

    :Params:

        - **name**: Name of the logger.
        - **logfile**, optional: Log file.
        - **console**, optional: Log to the console.
        - **maxlogsize**, optional: Maximal size of log file before rotating it.
        - **maxbackup**, optional: Maximal number of rotated files.
        - **sfmt**, optional: Format of log messages in log file.
        - **cfmt**, optional: Format of log message in console.
        - **asctime**, optional: Time format.
        - **level**, optional: Initialize logging level (see :meth:`set_loglevel`).
        - **colors**, optional: Use colors when formatting terminal messages?
        - **full_line**, optional: Colorize full line or just level name?
        - **redirect_warnings**, optional: Redirect messages issued by :mod:`warnings.warn`.
        - **redirect_stdout**, optional: Redirect messages issued to sys.stdout.
        - **redirect_stderr**, optional: Redirect messages issued to sys.stderr.

    :See also: :mod:`logging` module
    """
    def __init__(self, name, logfile=None, console=True, maxlogsize=0, maxbackup=0,
            cfmt='%(name)s [%(levelname)-8s] %(message)s',
            ffmt='%(asctime)s: %(name)s [%(levelname)-8s] %(message)s',
            asctime='%Y-%m-%d %H:%M',
            level='debug', colors=True, full_line=False,
            redirect_warnings=False, redirect_stdout=False, redirect_stderr=False):

        # Create or get logger
        self.logger = logger = logging.getLogger(name)

        # Handlers
        handlers = self.logger.handlers
        # - file
        if logfile is not None and logfile != '' and not any(
                [os.path.samefile(logfile,  l.baseFilename) for l in handlers
                    if isinstance(l, logging.handlers.RotatingFileHandler)]):
            checkdir(logfile, asfile=True)
            file =  logging.handlers.RotatingFileHandler(logfile,
                maxBytes=maxlogsize*1000, backupCount=maxbackup)
            file.setFormatter(logging.Formatter(ffmt, asctime))
            logger.addHandler(file)
        # - console
        if console and not any([(isinstance(l, logging.StreamHandler) and
                    not isinstance(l, logging.FileHandler)) for l in handlers]):
            console = logging.StreamHandler()
            if colors:
                console.setFormatter(ColoredFormatter(cfmt, full_line=full_line))
            else:
                console.setFormatter(logging.Formatter(cfmt))
            logger.addHandler(console)
        self.set_loglevel(level)

        # Redirections
        if redirect_warnings:
            warnings.showwarning = self.showwarning
        if redirect_stdout:
            if not isinstance(redirect_stdout, str):
                redirect_stdout = 'debug'
            sys.stdout = _Redirector_(getattr(self, redirect_stdout),
                prefix='STDOUT: ')
        if redirect_stderr:
            if not isinstance(redirect_stderr, str):
                redirect_stderr = 'warning'
            sys.stderr = _Redirector_(getattr(self, redirect_stderr),
                prefix='STDERR: ')


        # Announce
        logger.debug('*** Start log session ***')

    def debug(self, text, *args, **kwargs):
        """Send a debug message"""
        self.logger.debug(text, *args, **kwargs)

    def info(self, text, *args, **kwargs):
        """Send a info message"""
        self.logger.info(text, *args, **kwargs)

    def warning(self, text, *args, **kwargs):
        """Send a warning message"""
        self.logger.warning(text, *args, **kwargs)

    def showwarning(self, message, category, filename, lineno,
            file=None):
        self.warning(
            'REDIRECTED: %s:%s: %s:%s',
            filename, lineno,
            category.__name__, message,
        )


    def _log_and_exit_(self, slevel, text, *args, **kwargs):
        """Log a message and exit"""
        getattr(self.logger, slevel)(text, *args, **kwargs)
        mode = kwargs.pop('mode', None)
        if mode=='exiterr':
            mode = sys.exit
        elif mode=='exit':
            mode = 1
        if isinstance(mode, Exception):
            raise mode(text)
        elif callable(mode):
            mode(text)
        elif isinstance(mode, int):
            sys.exit(mode)


    def error(self, text, *args, **kwargs):
        """Send an error message"""
        self._log_and_exit_('error', text, *args, **kwargs)

    def critical(self, text, *args, **kwargs):
        """Send a critical message"""
        self._log_and_exit_('critical', text, *args, **kwargs)

    def set_loglevel(self, level=None, console=None, file=None):
        """Set the log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

        :Example:

            >>> logger.set_loglevel('DEBUG', console='INFO')
        """
        if level is not None:
            self.logger.setLevel(self._get_loglevel_(level))
        for handler in self.logger.handlers:
            if isinstance(handler, logging.handlers.RotatingFileHandler):
                if file is not None: handler.setLevel(self._get_loglevel_(file))
            elif console is not None and \
                isinstance(handler, logging.StreamHandler):
                handler.setLevel(self._get_loglevel_(console))

    def get_loglevel(self, asstring=False):
        """Get the log level as an integer or a string"""
        if asstring:
            for label in 'NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'FATAL':
                if self.logger.level == getattr(logging, label):
                    return label
            return 'NOTSET'
        return self.logger.level

    def _get_loglevel_(self, level):
        if level is None: level = 'debug'
        if isinstance(level, str):
            level = getattr(logging, level.upper(), 'DEBUG')
        return level

class TermColors(object):
    RED = ERROR = FAILURE = CRITICAL = '\033[31m'
    GREEN = INFO = OK = SUCCESS = '\033[32m'
    YELLOW = WARNING = '\033[33m'
    NORMAL = DEBUG = '\033[39m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    RESET = '\033[0m'

    def __init__(self):
        try: import curses
        except ImportError: curses = None
        import sys
        if curses:
            # This is a feature test which may end with an exception (eg. exec through ssh session, unset TERM env var)
            # We don't want to see this kind error (but we are masking other potential errors...)
            try: curses.setupterm()
            except: pass
        if not sys.stdout.isatty() or not curses or (curses.tigetstr('setf') is None and curses.tigetstr('setaf') is None):
            self.disable()

    def disable(self):
        for att in "RED ERROR FAILURE CRITICAL \
        GREEN INFO OK SUCCESS YELLOW WARNING \
        NORMAL DEBUG BLUE MAGENTA CYAN RESET".split():
            setattr(self, att, '')

    def format(self, text, color='NORMAL'):
        """Format a string for its color printing in a terminal

        - **text**: simple string message
        - *color*: color or debug level
        """
        color = color.upper()
        if not hasattr(self, color): return text
        return getattr(self, color)+text+self.RESET

def netcdf3():
    """Turn netcdf4 writing off with :mod:`cdms2`"""
    try:
        netcdf4(0, 0, 0)
    except:
        pass


def netcdf4(level=3, deflate=1, shuffle=1):
    """Turn netcdf4 writing on and suppress compression warning"""
    cdms2.setCompressionWarnings(0)
    cdms2.setNetcdfDeflateFlag(deflate)
    cdms2.setNetcdfShuffleFlag(shuffle)
    cdms2.setNetcdfDeflateLevelFlag(level)
netcdf4()

######################################################################
######################################################################

from .atime import ch_units, round_date, are_same_units, now, has_time_pattern, \
    tsel2slice, is_time, time_selector, itv_union, add_margin, strptime, \
    datetime as adatetime, comptime,  filter_time_selector, reltime
from .axes import get_checker, istime, islon, islat, islevel
from .color import land, simple_colors, get_cmap
from .grid import create_grid, get_xy, curv2rect, isgrid
from .grid.masking import polygons, convex_hull, rsamp, polygon_mask, create_polygon, \
    clip_shape
from .grid.regridding import griddata, xy2xy
from .grid.basemap import get_proj
from .misc import is_iterable,  broadcast, kwfilter, set_atts, create_selector, \
    squeeze_variable, checkdir, match_atts, get_atts, set_atts
from .phys.units import deg2m, m2deg
from .plot import map2, _colorbar_, savefigs as Savefigs, markers as Markers

