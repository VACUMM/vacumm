# -*- coding: utf8 -*-

"""
Misc tools

.. note::
    You can import it directly for example like this::

    >>> from vacumm.misc import auto_scale
    
"""
# Copyright or © or Copr. Actimar (contributor(s) : Stephane Raynaud) (2010)
# 
# raynaud@actimar.fr
# 
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

import numpy as N,MV2,cdtime, cdms2
MV = MV2
from cdms2 import createAxis, isVariable
from MV2 import nomask
import re,os
from copy import copy, deepcopy
from genutil import minmax
from types import *
MA=N.ma
from itertools import cycle
from genutil import grower
from vacumm import VACUMMError

__all__ = ['ismasked', 'bound_ops', 'auto_scale', 'basic_auto_scale', 'geo_scale', 
    'get_atts', 'cp_atts', 'set_atts', 'check_def_atts', 'iterable', 'isnumber', 
    'rm_html_tags', 'deg2str', 'lonlab', 'latlab', 'deplab', 'deg_from_dec', 'kwfilter', 
    'dict_filter','dict_aliases', 'dict_merge', 'mask_nan', 'write_ascii_time1d', 'xls_style',
    'FileTree', 'geodir', 'main_geodir', 'intersect', 'Att', 'broadcast', 'makeiter', 
    'get_svn_revision', 'dirsize', 'Cfg2Att', 'closeto', 'cp_props', 
    'zoombox','scalebox','history', 'dict_check_defaults', 'is_iterable', 
    'grow_variables', 'grow_depth', 'grow_lat',
    'create_selector', 'selector2str', 'split_selector', 'squeeze_variable', 'dict_copy_items', 
    "N_choose", 'MV2_concatenate', 'MV2_axisConcatenate', 'ArgList',
    'set_lang','set_lang_fr', 'lunique', 'tunique']
__all__.sort()

def broadcast(set, n, mode='last', **kwargs):
    """Broadcast ``set`` to the specified length ``n``
    
    :Params:
    
        - **set**: A single element or a sequence.
        - **n**: Final requested length.
        - **mode**, optional: Filling mode.
        
            - ``"last"``: Use the last element.
            - ``"first"``: Use the first element.
            - ``"cycle"``: Cycle through set.
            - if ``fillvalue`` is passed as a keyword, it is used to fill.
    
    :Example:
    
        >>> broadcast([2,4], 5)
        >>> broadcast(5, 4)
        >>> broadcast((2,3), 3)
        >>> broadcast((2,3), 3, mode='first')
        >>> broadcast([2,3], 3, mode='value', fillvalue=999)
    
    """
    set = makeiter(set)
    if n<=len(set): return set[:n]
    res = list(set)
    if kwargs.has_key('fillvalue'):
        fillvalue = makeiter(kwargs['fillvalue'])
    elif mode=='cycle':
        fillvalue = set
    elif mode=='first':
        fillvalue = set[:1]
    else:
        fillvalue = set[-1:]
    filliter = cycle(fillvalue)
    for i in xrange(n-len(set)):
        res.append(filliter.next())
    if isinstance(set, list): return res
    return res.__class__(res)

def makeiter(var):
    """Make var iterable as a list if not ietrable or a string"""
    if isinstance(var, str) or not hasattr(var, '__len__') or not callable(var.__len__) or \
        not hasattr(var, '__getslice__'):
        var = [var]
    return var

class Att(dict):
    """Class to create a dictionnary and access and set keys as attributes.
    
    You initialize and manage it as using :class:`dict`.
    
    :Example:
    
        >>> dd = Att(toto=3)
        >>> print dd['toto']
        3
        >>> print dd.toto
        3
        >>> dd.toto = 5
        >>> print dd['toto']
        5
    """
    
    def __getattr__(self, att):
        if self.has_key(att):
            return self[att]
        return dict.__getattribute__(self, att)

    def __setattr__(self, att, val):
        self[att] = val
        
def Cfg2Att(cfg):
    """Convert a :class:`~configobj.ConfigObj` object to an arborescence of :class:`~vacumm.misc.misc.Att` objects"""
    from configobj import ConfigObj, Section
    if isinstance(cfg, Att): return cfg
    assert isinstance(cfg, (Section, dict)), 'You must pass a ConfigObj object'
    a = Att()
    a.update(cfg.items())
    for sec in cfg.keys():
        if isinstance(cfg[sec], dict):
            a[sec] = Cfg2Att(cfg[sec])
    return a
    


def ismasked(arr):
    try:
        return ((arr.mask is None) or (arr.mask is nomask))
    except:
        return arr.mask is MA.nomask


def bound_ops(bounds):
    """Get operators that must be used for checking inclusion within boundaries.
    Returned operators (ops) must be used in the following way to return True if value is inside bounds:
    
    :Params:
    
        - **bounds**: Boundary closing indicator (like 'co')
    
    :Example:
    
        >>> ops[0](value,lower_bound)
        >>> ops[1](value,upper_bound)

    :Return: 
    
        A 2-element tuple of operators taken within (operator.ge,operator.gt,operator.le,operator.lt).
    """
    if bounds[0] == 'c':
        ops = (operator.ge, )
    else:
        ops = (operator.gt, )
    if bounds[1] == 'c':
        ops += (operator.le, )
    else:
        ops += (operator.lt, )
    return ops


def auto_scale(data = None, nmax = None,vmin = None, vmax = None,
               separators = None, fractions = False, symetric=False,
               keepminmax=False, **kwargs):
    """Computes levels according to a dataset and its range of values. Locators are on a 10-base. Different scaling can be used with this version.

    :Params:
    
    - **data**: The dataset
    - **vmax**, optional: Replaces max(data)
    - **vmin**, optional:   //   min(data)
    - **symetric**, optional: min_value and max_value are made symetric with respect to zero [default: False]
    - **nmax**, optional:  Maximal number of levels
    - **separators**, optional: Subdivision for irregular levels
    - **fractions**, optional: Consider separators as fractions (]0.,1.[)
    - **steps**, optional: Base 10 steps for finding levels [default: [1,2,2.5,5,10]]
    - **geo**, optional: Treat levels as geographical degrees [default: False]
    - Other parameters are given to MaxNLocator

    :Return: 
    
        Array of levels
    
    :See also:
    
        :func:`basic_auto_scale` :func:`geo_scale`
    """

    vmin = kwargs.get('min_value', vmin)
    vmax = kwargs.get('max_value', vmax)

    assert data is not None or (vmax is not None and vmin is not None), \
        'If data is not given, you must pass vmin and vmax as keywords'
        
    if data is not None:
        try: minv,maxv = minmax(data)
        except: minv,maxv = 0,1
    if vmin is not None:
        minv = vmin
    if vmax is not None:
        maxv = vmax

    if symetric:
        maxv = max([abs(maxv),abs(minv)])
        minv = -maxv
        
    if nmax is None: nmax = 7

    if separators is None:
        levels = basic_auto_scale(minv,maxv,nmax,**kwargs)

    else:

        try:
            separators.sort()
        except:
            separators = [separators,]
        nsep = len(separators)

        if fractions:
            relsep = []
            for s in separators:
                if s > 0. or s < 1.:
                    relsep.append(minv+(maxv-minv) * s)
            separators = relsep

        nint = nsep + 1
        nmax_sep = int(float(nmax)/float(nint))
        separators.insert(0,minv)
        separators.append(maxv)

        for i in xrange(nint):
            these_levels = basic_auto_scale(separators[i],
                                       separators[i+1],
                                       nmax_sep,**kwargs)
            if i == 0:
                levels = list(these_levels)
            else:
                last_level = levels[-1]
                for l in these_levels:
                    if l > last_level:
                        levels.append(l)
                        
    # Min/max
    if keepminmax:
        levels[0] = minv
        levels[-1] = maxv
        
    return levels

def basic_auto_scale(vmin,vmax,nmax=7,steps=[1,2,2.5,5,10],geo=False,minutes=False,**kwargs):
    """Computes levels according to a dataset and its range of values. Locators are on a 10-base.

    :Params:
    
        - **vmin/vmax**: Find levels around this range.
        - **nmax**, optional: Maximal number of steps.
        - **steps**, optional: Base 10 steps for finding levels [default: [1,2,2.5,5,10]]
        - **geo**, optional: Assume longitude or latitude degrees [default: False].
        - **minutes**, optional: If geo, find suitable levels to match nice minutes locations when locations have floating values (like 1.2) [default: False].
    
    :See also:
    
        :func:`auto_scale` :func:`geo_scale`
    """

    from matplotlib.ticker import MaxNLocator

    if geo:
        if (vmax-vmin) > 50.:
            steps = [1,2,3,6, 10]
        elif (vmax-vmin) > 3.:
            steps = [1,2,2.5,3,5,10]
        elif minutes:
#           # Minimal range
#           if (vmax-vmin) < 1:
#               vmax -= .5/60.
#               vmin += .5/60.
            # Minimal number of steps
            nmn = int(N.ceil((vmax-vmin)*60.))+1
#           nmax = min(max(nmax, int(N.ceil((vmax-vmin)*60.))+1), 10)
            nmax = 10
            # Steps
            if nmn > nmax:
                minutes_steps = [1, 10/6.,15/6.,20/6.,30/6.,10]
            else: # Strictly every minutes or less
                # TODO: must again use basic_auto_scale for minutes and seconds
                vmin = N.floor(vmin*60)
                vmax = N.ceil(vmax*60)
                if (vmax-vmin)<=2:
                    step = .5
                else:
                    step = 1.
                return N.arange(vmin, vmax+step, step)/60.
            # Try it
            myloc = basic_auto_scale(vmin,vmax,nmax,steps=minutes_steps,geo=False,**kwargs)
            # Check that we really have minutes, ie not only degrees
            if not N.allclose(N.array(myloc)%1., 0., atol=1.e-4):
                return myloc

    intv = (vmin, vmax)
    if nmax < 2: nmax = 2
    myloc = MaxNLocator(nmax+1,steps=steps,**kwargs)
    myloc.create_dummy_axis()
    myloc.set_view_interval(vmin, vmax)
    myloc.set_data_interval(vmin, vmax)

    locs = myloc()
    if len(locs)>1 and N.allclose(locs[1],  vmin): locs = locs[1:]
    if len(locs)>1 and N.allclose(locs[-2], vmax): locs = locs[:-1]
    return locs

def geo_scale(*args,**kwargs):
    """ :func:`auto_scale()` with geo=True 
    
    
    :See also:
    
        :func:`basic_auto_scale` :func:`auto_scale`
    """
    return auto_scale(geo=True,*args,**kwargs)


def get_atts(var, id=True, **kwargs):
    """ Get all attributes from a variable

    :Params:
    
        - **var**: Target variable
        - **id**, optional: Get also id [default: False]
        - Other keywords are set as attributes
    
    :See also:
    
        :func:`cp_atts` :func:`set_atts` :func:`check_def_atts`
    """
    atts = {}
    if hasattr(var,'attributes'):
        atts.update(var.attributes)
    if id and hasattr(var,'id'):
        atts['id'] = var.id
    atts.update(kwargs)
    return atts

def cp_atts(var1, var2, overwrite=True, select=None, exclude=None, **kwargs):
    """ Copy all atributes from one variable to another

    :Params:
    
        - **var1/var2**: Copy var1 attributes to var2
        - **id**, optional: Also copy id [default: False]
        - **overwrite**, optional: Overwrite attributes of ``var2``?
        - **select**, optional: Copy only these attributes. 
        - **exclude**, optional: Attributes that must be excuded from copy. 
        - Other keywords are set as attributes
    
    :See also:
    
        :func:`get_atts` :func:`set_atts`  :func:`check_def_atts`
    """
    # Default list
    atts = get_atts(var1,**kwargs)
    
    # Selection
    if atts and select is not None:
        if isinstance(select, basestring):
            select = [select]
        atts = dict(item for item in atts.items() if item[0] in select)
    
    # Exclusion
    if atts and exclude is not None:
        if isinstance(exclude, basestring):
            exclude = [exclude]
        atts = dict(item for item in atts.items() if item[0] not in exclude)
       
    # Set
    set_atts(var2, atts, overwrite=overwrite)

def set_atts(var, atts=None, overwrite=True, **kwargs):
    """Set attributes
    
    :Params:
    
        - **var**: Change attributes of var
        - **atts**, optional: A dictionary of attributes
        - **overwrite**, optional: Overwrite attributes of ``var2``?
        - Other keywords are set as attributes
    
    :See also:
    
        :func:`cp_atts` :func:`get_atts` :func:`check_def_atts`
    """
    if atts is None: atts = {}
    atts.update(kwargs)
    for att,val in atts.items():
        if overwrite or not hasattr(var, att):
            setattr(var, att, val)



def check_def_atts(obj, **defaults):
    """Check defaults attributes and set them if empty
    
    :See also:
    
        :func:`get_atts` :func:`cp_atts` :func:`set_atts` 
    """
    for att,val in defaults.items():
        if not hasattr(obj,att):
            setattr(obj,att,val)

def _pospos_(i, n):
    if isinstance(i, int):
        return i if i>=0 else (n-i)
    jj = []
    for j in i:
        jj.append(j if j>=0 else (n-j))
    return type(i)(jj)
    
def _negpos_(i, n):
    if isinstance(i, int):
        return i if i<0 else (i-n)
    jj = []
    for j in i:
        jj.append(j if j<0 else (j-n))
    return type(i)(jj)
    

def cp_props(var1, var2, axes=None, grid=True, atts=None, exaxes=None, exatts=None, owatts=True):
    """Copy properties of a variable to another variabes
    
    Proporties are attributes, axes and grid.
    
    :Params:
    
        - **var1/var2**: cdms variables.
        - **axes**, optional: Position of axes to copy. Positions are
          converted to relative position to that last dimension to handle
          the case of variables that don't have the same number of dimensions.
          Set it to ``False`` to not copy axes.
        - **exaxes**, optional: Position of axes to not copy.
        - **grid**, optional: Also copy the grid?
        - **atts**, optional: Attributes to copy. Set it to ``False`` to not
          copy attributes.
        - **exatts**, optional: Attributes to exclude from copy.
        - **owatts**, optional: Overwrite attributes?
        
    :Examples:
    
        >>> cp_props(var1, var2, atts=False)
        >>> cp_props(var1, var2, axes=[1,-1], atts=['units','long_name'])
        >>> cp_props(var1, var2, exaxes=0, grid=False)
        >>> cp_props(var1, var2, axes=False, exatts='units', owatts=False)
    """
    
    # Axes
    if axes is not False:
        if axes is None:
            axes = range(var1.ndim)
        axes = _negpos_(axes, var1.ndim)
        if exaxes is not None:
            if isinstance(exaxes, int):
                exaxes = [exaxes]
            exaxes = _negpos_(exaxes, var1.ndim)
            axes = [axis for axis in axes if axis not in exaxes]
        for iaxis in axes:
            if var1.ndim+iaxis <0 or var2.ndim+iaxis: continue
            var2.setAxis(iaxis, var1.getAxis(iaxis))
    
    # Grid
    if grid:
        grid = var1.getGrid()
        if grid is not None:
            from vacumm.misc.grid import set_grid
            set_grid(var2, grid)
            
    # Attributes
    if atts is not False:
        cp_atts(var1, var2, select=atts, exclude=exatts, overwrite=owatts)
        
    return var2

def is_iterable(obj, nostr=True, nogen=True):
    """Check if an object is iterable or not.

    :Params:
    
        - **obj**: Object to check

    :Return: True/False
    """
    
    if not nogen and type(obj) == types.GeneratorType: return True
    #try: len(obj)
    #except: return False
    if not (hasattr(obj, '__len__') and callable(obj.__len__)): return False
    if nostr: return not isinstance(obj, basestring)
    return True

iterable = is_iterable

def isnumber(var):
    return type(var) in [IntType,FloatType,LongType,ComplexType]



def rm_html_tags(str):
    """Remove html tags from a string

    :Params:
    
        - **str**: A string
    
    :Return: 
    
        Cleaned string
    
    Example:
    
        >>> rm_html_tags('<title>My title</title>')
        My title
    """

    import re
    return re.sub('<[^>]+>','',str)





_geolab_params = """- **decimal**, optional: Use decimal representation [default: True]
        - **fmt**, optional: Format for decimal representation [default: '%4.1f']
        - **no_seconds**, optional: Do not add seconds to degrees representation [default: False]
        - **tex**, optional: Use Tex degree symbol [default: False]
        - **no_symbol**, optional: Don't use degree symbol [default: False]
        - **no_zeros**, optional: Don't insert zeros when integer < 10 [default: False]
        - **auto**, optional: If True, find the ticks according to the range of values [default: False]
        - **auto_minutes**, optional: Automatically suppress degrees if value is not exactly a degree (just display minutes and seconds) else display degree [default: False]
    """

def deg2str(*args,**kwargs):
    return _geolab(*args,**kwargs)

def _geolab(vals,fmt='%.5g',longitude=True,decimal=True,tex=False,auto_minutes=False,
            no_seconds=False,no_symbol=False,no_zeros=False,auto=False, **kwargs):
    """Return a nice label for longitude or latitude

    Inspired from Basemap toolkit of Matplotlib
    """
    if kwargs.has_key('nosec'): no_seconds = kwargs['nosec']
    if kwargs.has_key('dec'): decimal = kwargs['dec']
    if kwargs.has_key('nosym'): no_symbol = kwargs['nosym']

    # Longitude/latitude
    if longitude:
        pstr,mstr = 'E','W'
    else:
        pstr,mstr = 'N','S'

    # Use tex strings
    if tex is None:
        from pylab import rcParams
        usetex = rcParams['text.usetex']
    elif tex:
        usetex = True
    else:
        usetex = False
    if usetex:
        sdeg = r"$^{\circ}$"
        smin = r"$^{'}$"
        ssec = r"$^{''}$"
    else:
        sdeg = u'\N{DEGREE SIGN}'
        smin = u"'"
        ssec = u"''"

    # Integer format
    if no_zeros:
        ifmt = '%i'
    else:
        ifmt = '%02i'
    gfmt = '%.2g'
    if decimal: auto_minutes = False

    # Loops on values
    if not iterable(vals):
        it = False
        vals = [vals,]
    else:
        it = True
    if auto:
        vals = minmax(vals)
    labs = []
    for val in vals:
        
        if val == 0.: # Equator or Greenwich
            if no_symbol:
                labstr = fmt
            else:
                labstr = fmt+sdeg
            labs.append(labstr%(val))

        else:
            if val > 180:
                val -= 360.
            elif val < -180.:
                val += 360.
            if val<0:
                sig = mstr
                val = -val
            else:
                sig = pstr
            nodeg = auto_minutes and abs(val%1) > 1.e-4
            if not decimal: # Minutes'(seconds'')
                fmt = ifmt
                if usetex and not ifmt.startswith('$'):
                    ifmt = r'$%s$'%ifmt
                dd,mm,ss = deg_from_dec(val)
                if no_seconds: # Round seconds to minute
                    mm = int(round(mm+ss/60.))
                    ss = 0
                if ifmt%ss != ifmt%0: # Seconds needed
                    vv = "%(ifmt)s%(smin)s%(gfmt)s%(ssec)s"% vars()
                elif ifmt%mm != ifmt%0: # Just minutes
                    vv = "%(ifmt)s%(smin)s" % vars()
                else: # Nothing
                    vv = ''
            else:
                vv=''
            if nodeg:
                labstr = vv
            elif no_symbol:
                labstr = fmt+vv
            else:
                labstr = fmt+sdeg+vv
            if not decimal:
                cc = []
                if not nodeg:
                    cc = [dd]
                if vv.endswith(ssec):
                    cc.extend([mm,ss])
                elif vv.endswith(smin):
                    cc.append(mm)
                labstr = labstr % tuple(cc)
            else:
                labstr = labstr % val
            if not nodeg: labstr += sig
#           if usetex:
#               labstr = r"%s${^'}" % val
#           else:
#               labstr = u'%s'%labstr
            labs.append(labstr)
    if it:
        return labs
    else:
        return labs[0]

def lonlab(longitudes,**kwargs):
    """Return nice longitude labels

    :Params:
    
        - **longitudes**: Value of longitudes   
        %s
    
    :See also:
    
        :func:`latlab` :func:`deplab` 
    """
    return _geolab(longitudes,**kwargs)
if lonlab.__doc__ is not None:
    lonlab.__doc__ = lonlab.__doc__ % _geolab_params

def latlab(latitudes,**kwargs):
    """Return nice latitude labels

    :Params:
    
        - **latitudes**: Value of latitudes
        %s
    
    :See also:
    
        :func:`lonlab` :func:`deplab` 
    """
    return _geolab(latitudes,longitude=False,**kwargs)
if latlab.__doc__ is not None:
    latlab.__doc__ = latlab.__doc__ % _geolab_params

def deplab(depths, fmt='%gm', auto=False, nosign=False):
    """Return well formatted depth labels

    :Params:
    
        - **depths**: Numerical depths

        - **fmt**, optional: Numeric format of the string (including units) [default: '%gm']
        - **auto**, optional: If True, find the ticks according to the range of depths [default: False]
        - **nosign**, optional: Absolute values are used.
    
    :See also:
    
        :func:`lonlab` :func:`latlab` 
    """

    single = not iterable(depths)
    if single: depths = [depths,]
    if auto:    depths = auto_scale(depths)

    labs = []
    for depth in depths:
        if nosign: depth = abs(depth)
        labs.append(fmt % depth)

    if single:
        return labs[0]
    else:
        return labs



def deg_from_dec(dec):
    """ Convert from decimal degrees to degrees,minutes,seconds """

    degrees = int(dec)
    dec = (dec-degrees)*60.
    minutes = int(dec)
    seconds = (dec-minutes)*60.
    
    # Approximations
    if abs(seconds) < 1.e-6:
        seconds = 0.
    elif (60.-seconds) < 1.e-6:
        seconds = 0.
        minutes += 1
        if minutes == 60:
            degrees += 1

    return degrees,minutes,seconds


def dict_filter(kwargs,filters, defaults=None, copy=False, short=False, keep=False, **kwadd):
    """Filter out kwargs (typically extra calling keywords)
    
    :Params:
    
        - **kwargs**: Dictionnary to filter.
        - **filters**: Single or list of prefixes.
        - *defaults*: dictionnary of default values for output fictionnary.
        - *copy*: Simply copy items, do not remove them from kwargs.
        - *short*: Allow prefixes to not end with ``"_"``.
        - *keep*: Keep prefix filter in output keys.
    
    :Example:
    
        >>> kwargs = {'basemap':'f', 'basemap_fillcontinents':True, 'quiet':False,'basemap_plot':False}
        >>> print kwfilter(kwargs,'basemap', defaults=dict(drawcoastlines=True,plot=True),good=True)
        {'plot': False, 'fillcontinents': True, 'good': True, 'basemap': 'f', 'drawcoastlines': True}
        >>> print kwargs
        {'quiet': False}
    """

    if isinstance(filters,str):
        filters = [filters]
    if copy:
        kwread = kwargs.get
    else:
        kwread = kwargs.pop

    # Set initial items
    kwout = {}
    for filter in filters:
        if not filter.endswith('_') and kwargs.has_key(filter):
            if isinstance(kwargs[filter],dict):
                kwout.update(kwread(filter))
            else:
                kwout[filter] = kwread(filter)
        if not short and not filter.endswith('_'): filter += '_'
        for att,val in kwargs.items():
            if att.startswith(filter) and att != filter:
                if keep:
                    kwout[att] = kwread(att)
                else:
                    kwout[att[len(filter):]] = kwread(att)

    # Add some items
    kwout.update(kwadd)

    # Set some default values
    if defaults is not None:
        for att,val in defaults.items():
            kwout.setdefault(att,val)
    return kwout

def kwfilter(*args, **kwargs):
    """Alias for :func:`dict_filter`"""
    return dict_filter(*args, **kwargs)

def dict_aliases(kwargs, aliases):
    """Remove duplicate entries in a dictionnary according to a list of aliases.
    The first alias has priority over the following.
    
    :Example:
    
        >>> kwargs = dict(min=4, title='Title', label='Label')
        >>> dict_aliases(kwargs, ['label', 'title'])
        {'min': 4, 'label': 'Label'}
    """
    if isinstance(aliases, (str, unicode)):
        aliases = [aliases]
    ref_alias = aliases[0]
    for i, alias in enumerate(aliases):
        if kwargs.has_key(alias):
            if i:
                kwargs[ref_alias] = kwargs.pop(alias)
            for j in xrange(i+1, len(aliases)):
                if kwargs.has_key(aliases[j]):
                    del kwargs[aliases[j]]
            break
    return kwargs

def dict_check_defaults(kwargs, **defs):
    """Check that a dictionary has some default values"""
    if defs is None: defs = {}
    for item in defs.iteritems():
        kwargs.setdefault(*item)
    return kwargs
    
def lunique(mylist):
    """Uniquify a list to a new list"""
    ulist = []
    seen = {}
    for item in mylist:
        if item in seen: continue
        seen[item] = True 
        ulist.append(item)
    return ulist
    
def tunique(mytuple):
    """Uniquify a tuple to a new tuple"""
    utuple = []
    seen = {}
    for item in mytuple:
        if item in seen: continue
        seen[item] = True 
        utuple += item, 
    return utuple
       

def dict_merge(*dd, **kwargs):
    """Merge dictionaries
    
    :Params:
    
        - **dd**: Argument ar interpreted as dictionary to merge.
          Those who are not dictionaries are skipped.
        - **mergesubdicts**, optional: Also merge dictionary items 
          (like in a tree) [default: True].
        - **mergetuples**, optional: Also merge tuple items [default: False].
        - **mergelists**, optional: Also merge list items [default: False].
        - **unique**, optional: Uniquify lists and tuples [default: True].
        - **skipnones**, optional: Skip Nones [default: True].
     
    """
    # Options
    mergesubdicts = kwargs.get('mergesubdicts', True)
    mergelists = kwargs.get('mergelists', False)
    mergetuples = kwargs.get('mergetuples', False)
    unique = kwargs.get('unique', True)
    skipnones = kwargs.get('skipnones', True)
   
    # Loop
    outd = {}
    for d in dd:
        if not isinstance(d, dict): continue
        for key, val in d.iteritems():
            if skipnones and val is None: continue
            if key not in outd: # Not set so we set
                outd[key] = val
            elif mergesubdicts and isinstance(outd[key], dict) and isinstance(val, dict): # Merge subdict
                outd[key] = dict_merge(outd[key] , val, **kwargs)
            elif mergelists and isinstance(outd[key], list) and isinstance(val, list): # Merge lists
                outd[key] += val
                if unique:
                    outd[key] = lunique(outd[key])
            elif mergetuples and isinstance(outd[key], tuple) and isinstance(val, tuple): # Merge tuples
                outd[key] += val
                if unique:
                    outd[key] = tunique(outd[key])
                
    return outd

def dict_copy_items(ddi, ddo, keys):
    """Copy existing items of an array to another one
    
    :Params:
    
        - **ddi**: Dictionary from which items are taken.
        - **ddo**: Dictionary in which items are put.
        - **keys**: single or list of keys.
    """
    if isinstance(keys, basestring):
        keys = [keys]
    single = isinstance(ddo, dict)
    if single: ddo = [ddo]
    for key in keys:
        if key in ddi:
            for dd in ddo:
                dd[key] = ddi[key]
    if single: return ddo[0]
    return ddo
    
def mask_nan(input):
    """Mask NaN from a variable
    
    :Params:
    
        - **input**: Array (numpy, numpy.ma or MV2)
    
    .. note::
    
        If input is pure numpy, it is converted to numpy.ma
    
    :Example:
    
        >>> var = N.array([1,N.nan])
        >>> print mask_nan(var)
        [1.0 --]
    """
    isma = isVariable(input) or MA.isMA(input)
    if isma:
        pure_num = input.filled()
    else:
        pure_num = N.asarray(input)
    output = MA.masked_where(N.isnan(pure_num), input, copy=False)
    if isma:
        input[:] = output
    return output
    
#   # Masking function
#   if isVariable(input):
#       M = MV
#   else:
#       M = MA
#
#   # Work on numerical values
#   if MA.isMA(input):
#       num = MA.filled(input)
#   else:
#       num = input
#
#   # Compute the mask
#   a=N.greater(num,0.)
#   b=N.less_equal(num,0.)
#   mask = N.logical_not(N.logical_or(a,b))
#   del a,b
#
#   # Apply the mask if needed
#   # Output is always at least a MA class
#   if N.sometrue(mask):
#       output = M.array(input, mask= mask)
#   else:
#       output = M.array(input, copy=0)
#
#   return output
    
def write_ascii_time1d(var,file,fmt='%g'):
    """Write an ascii file in the following format: YYY/MM/DD HH:MN:SS DATA where DATA is one column.
    
    :Params:
    
        - **var**: a cdms variable WITH A TIME AXIS
        - **file**: output file name
        - **fmt**, optional: format of DATA [default: '%g']
    """

    import vacumm.misc.atime as T
    try:
        time = var.getTime().asComponentTime()
    except:
        raise Exception,'[write_ascii1d] No axis time attached to this variable'

    line_fmt = '%s '+fmt

    f = open(file,'w')
    var.setMissing(999.)
    var = MV.filled(var)
    for it in xrange(len(time)):
        stime = T.num_to_ascii(time[it])
        f.write(line_fmt % (stime,var[it].flat[0])+'\n')
    f.close()

def xls_style(style=None,b=None,i=None,u=None,c=None,o=None,bd=None,fmt=None,va=None,ha=None, f=None, s=None, n=None, copy=True, **kwargs):
    """Excel style sheet for pyExcelerator cell objects
    
    :Params:

        - **style**, optional: Style object to update.
        - **scopy**, optional: Make a copy of ``style`` if passed.
        - **b**, optional: Bold [True/False].
        - **i**, optional: italic[True/False].
        - **u**, optional: Underline [True/False].
        - **c**, optional: Colour index [from 0 to 82!, with 0=black, 1=white, 2=red, 3=green, 4=blue, 5=yellow, 6=magenta, 7=cyan].
        - **o**, optional: Outline [True/False].
        - **bd**, optional: Borders [True/int/False/dict(top=int...)].
        - **fmt**, optional: Format ['general','0.0',...].
        - **va**, optional: Vertical alignment ['top','bottom','center','justified'].
        - **ha**, optional: Horizontal alignment ['left','right','center','justified'].
        
    :Example:
    
        >>> style_num = xls_style(fmt='0.00')
        >>> style_num_left = xls_style(style_num, left=2, copy=True)
        >>> style_full = xls_style(c=2, b=True, va='center', ha='justified')
        
    :Tutorial: :ref:`user.tut.misc.io.xls`
    """
    from pyExcelerator import XFStyle,Borders,Alignment, Font
    # Style object
    if style is None:
        style = XFStyle()
    elif copy:
        style = deepcopy(style)
    # Format
    if fmt is not None:
        style.num_format_str = fmt
    # Font syle
    if not hasattr(style, 'font'):
        style.font = Font()
    if b is not None:
        style.font.bold = bool(b)
    if i is not None:
        style.font.italic = bool(i)
    if u is not None:
        style.font.underline = bool(u)
    if o is not None:
        style.font.outline = bool(o)
    if c is not None:
        style.font.colour_index = c#hex(c)
    if n is not None:
        style.font.name = n
    if s is not None:
        style.font.height = int(s)
    # Borders
    bdnames = 'top', 'left', 'right', 'bottom'
    if bd is not None and not isinstance(bd, dict):
        bd = int(bd)
        bd = dict(bottom=bd,top=bd,left=bd,right=bd)
    else:
        bd = {}
    for bdname in bdnames:
        if kwargs.has_key(bdname):
            value = kwargs.pop(bdname)
            if value is not None:
                bd[bdname] = int(value)
    if not hasattr(style, 'borders'):
        style.borders = Borders()
    for ik, (key,val) in enumerate(bd.items()):
        setattr(style.borders,key,val)
    # Alignment
    if not hasattr(style, 'alignment'):
        style.alignment = Alignment()
    if va is not None:
        style.alignment.vert = getattr(Alignment,'VERT_'+va.upper())
    if ha is not None:
        style.alignment.horz = getattr(Alignment,'HORZ_'+ha.upper())
    return style





class FileTree(object):
    """Build a file tree

    :Params:
    
        - **input_dir**: Input directory.
        - **patterns**, optional: A string (or a list of strings) indicating which REGULAR EXPRESSION patterns files must match (using glob.glob) [default: '.*']
        - **exclude**, optional:: A string (or a list of strings) indicating which REGULAR EXPRESSION patterns files must not match (using re.search) [default: ['CVS/','.svn/']]
        - **include**, optional:: A string (or a list of strings) indicating which REGULAR EXPRESSION patterns files must match if they are excluded by 'exclude' patterns (using re.search) [default: None]
        - **relative**, optional:: Return file names relative to input_dir [default: False]
    """
    default_patterns = dict(patterns=['.*'],exclude=['~$','^CVS$/','^\.svn$/','^\.DS_Store$','\.db$','\.ini$'],include=[],)

    def __init__(self,input_dir,relative=False,scan=True,**kwargs):

        if not os.path.exists(input_dir):
            raise 'Directory not found: '+input_dir
        self._input_dir = input_dir
        self._relative = relative

        self._set_selectors_(patterns=True,exclude=True,include=True,update=False)
        self._tree = None
        if scan: self.scan(**kwargs)


    def file_list(self,**kwargs):
        if self._tree is None:
            self.scan(**kwargs)
        return self._tree

    def __call__(self,**kwargs):
        return self.file_list(**kwargs)

    def __str__(self,sep='\n'):
        return sep.join(self._tree)

    def scan(self,**kwargs):
        """Recursive scan to list files"""
        kwargs['update'] = False
        self._set_selectors_(**kwargs)
        self._tree = []
        for current_dir,dirs,files in os.walk(self._input_dir):
            for isfile,targets in enumerate((copy(dirs),copy(files))):
                for target in targets:
                    if (isfile and not self._check_target_(current_dir,'patterns',target)) or \
                        (self._check_target_(current_dir,'exclude',target) and \
                         not self._check_target_(current_dir,'include',target)) :
                        if not isfile:
                            dirs.remove(target)
                    else:
                        if isfile:
                            this_file = os.path.join(current_dir,target)
                            if self._relative:
                                this_file = this_file[len(self._input_dir)+1:]
                            self._tree.append(this_file)

    def _check_target_(self,current_dir,pattern_type,target):
        for pat in getattr(self,'_'+pattern_type):
            thistarget = target
            if pat.pattern.endswith('/$') :
                if not os.path.isdir(os.path.join(current_dir,target)):
                    # Target must be a directory
                    continue
                else:
                    thistarget += '/'
            if pat.pattern.startswith('^'+self._input_dir):
                # Pattern is matched against full path to root dir
                thistarget = os.path.join(current_dir,thistarget)
            if pat.search(thistarget) is not None:
                return True
        return False


    def set_patterns(self,value,**kwargs):
        self._set_selector_('patterns',value,**kwargs)
    def set_exclude(self,value,**kwargs):
        self._set_selector_('exclude',value,**kwargs)
    def set_include(self,value,**kwargs):
        self._set_selector_('include',value,**kwargs)

    def _set_selector_(self,seltype,values,update=True,append=False):
        if values is None:
            values = []
        elif isinstance(values,str):
            values = [values]
        else:
            values = list(values)
        if len(values) and values[0] == '+':
            append = True
            del values[0]
        if append:
            vals = setattr(self,'_'+seltype)
        else:
            vals = []
        for val in values:
            if val.startswith('//'):
                val = '^'+os.path.join(self._input_dir,val[2:])
            if val.endswith('/'): val += '$'
            vals.append(re.compile(val))
        setattr(self,'_'+seltype,vals)
        if update: self.scan()

    def _set_selectors_(self,update=True,**kwargs):
        for key,val in kwargs.items():
            if key not in self.default_patterns.keys(): continue
            if val is True: val = self.default_patterns[key]
            self._set_selector_(key,val,False)
        if update: self.scan()

def geodir(direction, from_north=True):
    """Return a direction in degrees in the form 'WNW'"""
    if from_north:
        direction = 90. - direction
    direction = (direction+360) % 360.
    labels = ['E','ENE','NE','NNE','N','NNW','NW','WNW','W','WSW','SW','SSW','S','SSE','SE','ESE']
    nlab = len(labels)
    dtheta = 360./nlab
    return labels[int((direction+dtheta/2)/dtheta)%nlab]


def main_geodir(directions, amp=None, num=False, res=22.5, getamp=False, **kwargs):
    """Return the dominant direction from a set of directions in degrees"""
    cartesian = isinstance(directions,tuple)
    if cartesian:
        kwargs['from_north'] = False
        directions = list(directions)
        for i in 0,1:
            directions[i] = N.ravel(directions[i])
        if amp is None:
            amp = N.sqrt(directions[0]**2+directions[1]**2)
    else:
        directions = N.ravel(N.array(directions).astype('f'))
    if not num and amp is None:
        assert int(360./res) == 360./res,'360. must be a multiple of your resolution: %g'%res
        idirections = (((directions+360.+res/2) % 360.)/res).astype('i').tolist()
        bincount = N.zeros(int(360./res))
        for id in idirections:
            bincount[id] += 1
        return geodir(res*N.argmax(bincount),**kwargs)
    else:
        if amp is None:
            amp = 1.
        elif not cartesian:
            amp = N.ravel(N.array(amp).astype('f'))
        else:
            xmean = N.average(directions[0])
            ymean = N.average(directions[1])
        if not cartesian:
            xmean = N.average(amp*N.cos(directions*N.pi/180.))
            ymean = N.average(amp*N.sin(directions*N.pi/180.))
        assert xmean != 0 or ymean != 0, 'There is no dominant directions'
        md = N.arctan2(ymean,xmean)*180./N.pi
        if not num:
            return geodir(md,**kwargs)
        if getamp:
            return md,N.sqrt(xmean**2+ymean**2)
        return N.arctan2(ymean,xmean)*180./N.pi


def intersect(seg1,seg2,length=False):
    """Intersection of two segments
    
    :Example:
        
        >>> intersect([1,10],(5,12))
        [5, 10]
    """
    if max(seg1) < min(seg2) or min(seg2) > max(seg1):
        if length: return 0.
        return None
    x0,x1 = max(min(seg1),min(seg2)),min(max(seg1),max(seg2))
    if length: return x1-x0
    return type(seg1)([x0,x1])

def get_svn_revision(path, max=False):
    """Get the revision number of a path
    
    Adapted from :class:`numpy.disutils.misc_util.Configuration`
    """
    revision = None
    if not max:
        import sys
        m = None
        try:
            sin, sout = os.popen4('svnversion')
            m = re.match(r'(?P<revision>\d+)', sout.read())
        except:
            pass
        if m:
            revision = int(m.group('revision'))
            return revision
        from numpy.distutils.misc_util import njoin
        if sys.platform=='win32' and os.environ.get('SVN_ASP_DOT_NET_HACK',None):
            entries = njoin(path,'_svn','entries')
        else:
            entries = njoin(path,'.svn','entries')
        if os.path.isfile(entries):
            f = open(entries)
            fstr = f.read()
            f.close()
            if fstr[:5] == '<?xml':  # pre 1.4
                m = re.search(r'revision="(?P<revision>\d+)"',fstr)
                if m:
                    revision = int(m.group('revision'))
                else:  # non-xml entries file --- check to be sure that
                    m = re.search(r'dir[\n\r]+(?P<revision>\d+)', fstr)
                    if m:
                        revision = int(m.group('revision'))
    # SR hack
    if revision is None:
        import subprocess
        path = os.path.abspath(path)
        revision = subprocess.Popen(["svnversion", path], stdout=subprocess.PIPE).communicate()[0].split(':')[-int(max)].split('\n')[0]
    return revision

def dirsize(folder, units='b'):
    """Get the size of a directory
    
    :Params:
    
        - **folder**: path of the directory
        - **units**, optional:
    
            - ``"b"``: bytes
            - ``"k"``: Kb
            - ``"m"``: Mb
            - ``"g"``: Gb
            - ``"t"``: Tb
        
    """
    assert os.path.exists(folder), 'Path not found: '+folder
    if not os.path.isdir(folder):
        path = os.path.dirname(folder)
    folder_size = 0
    for (path, dirs, files) in os.walk(folder):
        for file in files:
            filename = os.path.join(path, file)
            folder_size += os.path.getsize(filename)
    units = str(units).lower()
    uu = 'kmgt'.find(units[0])
    return folder_size/(1024.**(1+uu))
    
def closeto(a, b, rtol=1.e-5, atol=1.e-8):
    """Check which values of a numeric array are close to those of another array
    
    :Params:
    
        - **a**: A :class:`numpy.ndarray` variable.
        - **b**: Another array.
        
    :See also: :func:`~numpy.allclose`
    """
    mask = N.ma.nomask
    if N.ma.isMA(a):
        mask = a.mask
        x = a.filled()
    else:
        x = N.array(a, copy=False)
    if N.ma.isMA(b):
        mask |= b.mask
        y = b.filled()
    else:
        y = N.array(b, copy=False)
    xinf = N.isinf(x)
    if not xinf.any():
        res = N.less_equal(N.absolute(x-y), atol + rtol * (N.absolute(y)))
    else:
        b = x==y
        g = ~xinf
        x = x[~xinf]
        y = y[~xinf]
        b[~xinf] = N.less_equal(N.absolute(x-y), atol + rtol * N.absolute(y))
        res = b
    if mask is not N.ma.nomask:
        res &= ~mask
    return res

    
def MV2_concatenate(arrays, axis=0, axisid=None, axisattributes=None, copy=True):
    if cdms2.isVariable(arrays):
        if copy: arrays = arrays.clone()
        return arrays
    if len(arrays)==1:
        if copy: arrays = [arrays[0].clone()]
        return arrays[0]
    var = MV2.concatenate(arrays, axis=axis, axisid=None, axisattributes=None)
    var.setAxis(0, MV2_axisConcatenate([v.getAxis(0) for v in arrays], 
        id=axisid, attributes=axisattributes, copy=copy))
    if arrays:
        cp_atts(arrays[0], var)
    return var
MV2_concatenate.__doc__ = MV2.concatenate.__doc__

def MV2_axisConcatenate(axes, id=None, attributes=None, copy=True):
    """Advanced version of MV2.axisConcatenate"""
    from axes import isaxis
    # Single
    if isaxis(axes):
        if copy: axes = axes.clone()
        return axes
    if len(axes)==1:
        if copy: axes = [axes[0].clone()]
        return axes[0]
    
    # Attributes
    if id is None: id = axes[0].id
    if attributes is None: attributes = axes[0].attributes
    
    # Time units
    if axes[0].isTime():
        units = getattr(axes[0], 'units', None)
        if units is not None:
            for axis in axes[1:]:
                un = getattr(axis, 'units', None)
                if un is not None and un!=units:
                    axis.toRelativeTime(units)
    
    # Basic concatenate
    axis = MV2.axisConcatenate(axes, id=id, attributes=attributes)
    return axis
MV2_axisConcatenate.__doc__ = MV2.axisConcatenate.__doc__


def scalebox(box, factor):
    """Alter box limits with a zoom factor
    
    :Params:
    
        - **box**: A list of ``[xmin,ymin,xmax,ymax]``.
        - **factor**: Zoom factor where 1 means no change,
          and a factor < 1 means a smaller box.
          
    :Example:
    
        >>> scalebox([0,0,1,2], 1.1)
        [-0.55, -1.10, 1.55, 2.55]
    """
    if cdms2.isVariable(box):
        if box.getGrid() is None: 
            raise TypeError('You must provide a variable with a valid grid')
        box = box.getGrid()
    from grid.misc import isgrid, get_xy
    if isgrid(box):
        lon, lat = get_xy(box, num=True)
        xmin = lon.min()
        xmax = lon.max()
        ymin = lat.min()
        ymax = lat.max()
    else:
        xmin, ymin, xmax, ymax = box
    dx = (xmax-xmin)*(factor-1)*.5
    dy = (ymax-ymin)*(factor-1)*.5
    newbox = xmin-dx, ymin-dy, xmax+dx, ymax+dy
    if isinstance(box,(list, tuple)):
        return box.__class__(newbox)
    return N.array(newbox)
    
def zoombox(box, factor):
    """Alias for :func:`scalebox` with ``1/factor`` as zoom factor."""
    return scalebox(box, 1/factor)

def history(nbcommand=None):
    """Display the command history for the interactive python."""
    import readline
    if nbcommand != None:
        index = range(readline.get_current_history_length()-nbcommand,readline.get_current_history_length())
    else:
        index = range(readline.get_current_history_length())
    
    for i in index:
        print readline.get_history_item(i)


def create_selector(*args, **kwargs):
    """Create a :class:`cdms2.selectors.Selector`
    
    :Params:
    
        - **args**: Each item of args is treated in a special way:
        
            - A dictionary is treated using: ``selector.refine(**args[i])``
            - A list  is treated using: ``selector.refine(*args[i])``
            - Else (single slice, tuple, etc): ``selector.refine(args[i])``
        
        - **kwargs**: Items are integrated using: ``selector.refine(**kwargs)``
        
    :Examples:
    
        >>> create_selector((time1,time2),slice(0,1),lat=(lat1,lat2))
        >>> create_selector([(time1,time2),(lat1,lat2)],(lon1,lon2))
        >>> selector2 = create_selector(selector1, lon=(lon1,lon2))
    """
    selector = cdms2.selectors.Selector()
    for arg in args:
        if arg is None: continue
        if isinstance(arg, dict):
            selector.refine(**arg)
        elif isinstance(arg, list):
            selector.refine(*arg)
        else:selector.refine(arg)
    if kwargs:
        for key, val in kwargs.items():
            if val is not None:
                selector.refine(**{key:val})
    return selector
        
def selector2str(selector):
    """Convert a selector to a generic string where <obj ...> labels are removed
    
    :Example:
    
        >>> sel = cdms2.selectors.Selector(slice(2,5),time=(4,6))
        >>> selector2str(sel)
        '((slice(2, 5, None)),("time", (4, 6)))'
    """
    if not isinstance(selector, cdms2.selectors.Selector):
        selector = create_selector(selector)
    ss = []
    for c in selector.components():
        s = str(c)
        s = s[s.index('>')+1:]
        ss.append(s)
    return '('+','.join(ss)+')'

def split_selector(selector):
    """Split a selector into (positionalComponents, axisComponents).
    
    :Params:
        - **selector**: cdms2.selectors.Selector
    
    :Return: a tuple of:
        - the list of positionalComponents
        - the dict of axisComponents
    
    :Exemple:
    
        >>> selector = cdms2.selectors.Selector(slice(2,5),time=(4,6))
        >>> split_selector(selector)
        ((slice(2, 5, None),), {'time': (4, 6)})
    
    """
    posed = []
    named = {}
    for c in selector.components():
        if isinstance(c, cdms2.selectors.positionalComponent):
            posed.append(c.v)
        else:
            named[c.id] = c.spec
    return tuple(posed), named

def squeeze_variable(var, spec=True):
    """Remove singleton axes from a MV2 variable
    
    :Params:
    
        - **var**: MV2 array.
        - **spec**, optional: Squeeze specification.
        
            - ``True`` or ``None``: Simply squeeze out all singletons.
            - ``False``: Does nothing (return the same variable).
            - A string containing geo letters like 'xyzt' to remove
              axes according to their type.
              
              
    :Examples:
    
        >>> squeeze_variable(var)
        >>> squeeze_variable(var, spec='tx')
    """
    # Nothing to do
    if spec is False or spec==0: return var
    
    # Squeeze all singletons
    if spec is True or spec is None or spec==1:
        return var(squeeze=1)
        
    # Squeeze selection
    if isinstance(spec, basestring):
        for axtype in spec:
            order = var.getOrder()
            if not axtype in order: continue
            iaxis = order.index(axtype)
            if len(var.getAxis(iaxis))==1:
                ss = [':']*var.ndim
                ss[iaxis] = 0
                var = var[tuple(ss)]
    return var

def grow_variables(var1, var2):
    """Grow dimensions of var1 and var2 until their match
    
    It is an improved version :func:`genutil.grower.grower`.
    """
    # Guess and merge orders
    from axes import merge_orders, set_order
    from grid.misc import set_grid
    order1, order2 = merge_orders(var1, var2)
    if '-' in order1+order2:
        raise VACUMMError("All axes must have a known type (xyzt) to grow variables")
    var1 = MV2.asarray(var1)
    var2 = MV2.asarray(var2)
    set_order(var1, order1)
    set_order(var2, order2)
    rev = len(order1)<len(order2)
    if rev: # reverse?
        var1, var2 = var2, var1
    
    # Grow variables
    grid1 = var1.getGrid()
    grid2 = var2.getGrid()
    var1, var2 = grower(var1, var2)
    if grid1 is not None or grid2 is not None:
        ggrid1 = grid1 or grid2
        ggrid2 = grid2 or grid1
        set_grid(var1, ggrid1)
        set_grid(var2, ggrid2)
    if rev: # reverse?
        var1, var2 = var2, var1
    return var1, var2

def grow_depth(var, depth=None, mode=None, default=0., getvar=True):
    """Make depth and variable have the same shape
    
    :Params:
        
        - **var**: A MV2 array with a depth dimension.
        - **depth**, optional: Depth axis or array. If None,
          it is guessed from var. If not found, it defaults to ``default``.
        - **default**, optional: Default value for depth if not found.
        - **mode**, optional: How to handle case where shape of var 
          is changed:
          
            - ``None``: Nothing happens
            - ``"warn"``: A warning with :func:`~warnings.warn`.
            - ``"raise"``: Raise a :exc:`~vacumm.VACUMMError` exception.
    
    :Return: ``var, depth`` where ``depth`` may be a scalar equal to ``default``
    """
    if depth is None:
        depth = var.getLevel()
        if depth is None:
            depth = default
    if not N.isscalar(depth):
        depth = MV2.asarray(depth)
        if depth.ndim==1:
            depth.getAxis(0).designateLevel()
        if depth.shape!=var.shape:
            s = var.shape
            var, depth = grow_variables(var, depth)
            if var.shape!=s:
                msg = 'Shape of variable has changed (%s->%s) when guessing depth'%(s, var.shape)
                if mode=='warn':
                    warn(msg)
                elif mode=='raise':
                    raise VACUMMError(msg)
    if getvar: return var, depth
    return depth
    
def grow_lat(var, lat=None, mode=None, default=None, getvar=True):
    """Make latitude and variable have the same shape
    
    :Params:
        
        - **var**: A MV2 array with a lat dimension.
        - **depth**, optional: Latitude axis or array. If None,
          it is guessed from var.
        - **default**, optional: Default value for lat if not found.
        - **mode**, optional: How to handle case where shape of var 
          is changed:
          
            - ``None``: Nothing happens
            - ``"warn"``: A warning with :func:`~warnings.warn`.
            - ``"raise"``: Raise a :exc:`~vacumm.VACUMMError` exception.
    
    :Return: ``var, lat where lat may be ``None``
    """
    if lat is None:
        lat = var.getLatitude()
        if lat is None:
            lat = default
    if lat is not None and not N.isscalar(lat):
        lat = MV2.asarray(lat)
        lat.getAxis(0).designateLatitude()
        if lat.ndim==2:
            lat.getAxis(1).designateLongitude()
        s = var.shape
        var, lat = grow_variables(var, lat)
        if var.shape!=s:
            msg = 'Shape of variable has changed (%s->%s) when guessing latitude'%(s, var.shape)
            if mode=='warn':
                warn(msg)
            elif mode=='raise':
                raise VACUMMError(msg)
    if getvar: return var, lat
    return lat


def N_choose(a, choices, out=None, mode='raise'):
    """Robust version to :func:`numpy.choose` valid also
    for masked arrays.
    
    The func:`numpy.choose` function fails when length of 
    choices is greater than 32!
    You can verify it with::
    
        >>> numpy.choose([0], numpy.zeros((50, 1))) ! FAIL
        >>> N_choose([0], numpy.zeros((50, 1))) ! SUCCESS
        
    .. note:: Performances of this function can be improved, even
        if it should not take alot of time to run in the vast
        majority of cases.
    """
    a = N.asarray(a)
    if mode=='wrap':
        a = a % len(choices)
    elif mode=='clip':
        a = N.clip(a, 0, len(choices)-1)
    if out is None:
        out = choices[0].copy()+100
    N_where = N.ma.where if N.ma.isMA(choices[0]) else N.where
    for i in xrange(a.max()+1):
        out[:] = N_where(a==i, choices[i], out)
    return out

class ArgList(object):
    """Utility to manage always arguments as list and return results as input
    
    :Examples:
    
        >>> a = 'a'
        >>> al = ArgList(a)
        >>> al.get() # input for function as list
        ['a']
        >>> al.put(['aa']) # output as input
        'aa'
        
        >>> a = ['a','b']
        >>> al = ArgList(a)
        >>> al.get()
        ['a', 'b']
        >>> al.put(['aa'])
        ['aa']
        
    """
    def __init__(self, argsi):
        self.single = not isinstance(argsi, list)
        self.argsi = argsi
    def get(self):
        return [argsi] if self.single else argsi
    def put(self, argso):
        so = not isinstance(argso, list)
        if (so and self.single) or (not so and not self.single):
            return argso
        if so and not self.single:
            return [argso]
        return argso[0]
 
def set_lang(default='en_US.UTF-8', num=None):
    """Set the default and numeric languages

    The numeric language defaults to the default language
    """
    import locale
    if num is None: num = base
    os.environ['LANG'] = default
    locale.setlocale(locale.LC_ALL, default)
    os.environ['LC_NUMERIC'] = num
    locale.setlocale(locale.LC_NUMERIC, num)
    
def set_lang_fr(ennum=True):
    """Set lang to french, except the numeric lang wich is set to en by default"""
    num = 'en_US.UTF-8' if ennum else 'fr_FR.UTF-8'
    set_lang(default='fr_FR.UTF-8', num=num)

#from phys.units import *
#from phys.constants import *



