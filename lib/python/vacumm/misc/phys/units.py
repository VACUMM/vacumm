# -*- coding: utf8 -*-
"Unit conversions"

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

import re
from constants import *
import numpy as N
import MV2

__all__ = [
    'kt2ms', 'ms2kt', 'deg2m', 'm2deg', 'ms2bf', 'dms2deg', 'deg2dms',
    'mph2ms', 'ms2mph', 'tometric', 'kel2degc', 'degc2kel',
    'convert_units', 'basic_proj', 'uuconvert',
    'rad2deg', 'deg2rad', 'vect2mod', 'vect2dir', 'vect2moddir',
    'moddir2vectx', 'moddir2vecty', 'moddir2vectxy',
    'strfsize', 'strpsize', 'uuconvert',
]

def uuconvert(value, oldunits, newunits):
    """Change units using unidata units connverter"""
    try:
        from genutil.udunits import udunits_wrap
    except:
        from unidata import udunits_wrap
    s, i = udunits_wrap.convert(oldunits, newunits)
    return value * s + i
convert_units = uuconvert

############################################################
def kt2ms(nd):
    """Convert nds to m/s"""
    return uuconvert(nd, 'kt', 'm/s')

def ms2kt(ms):
    """Convert m/s to nds"""
    return uuconvert(ms, 'm/s', 'kt')

ms2nd = ms2kt
nd2ms = kt2ms

def mph2ms(mph):
    """Convert from mph to m/s"""
    return 0.44704*mph

def ms2mph(ms):
    """Convert from m/s to mph"""
    return 2.237*ms

############################################################
def deg2m(degrees, lat=None):
    """Convert a distance in degrees to meters

    - **degrees**: Distance in degrees
    - *lat*: optional latitude, defaults to 0.

    Return: Distance in meters

    .. seealso::

        :func:`m2deg`
    """

    if lat is None:
        lat = 0.
    else:
        lat = N.array(lat)

    return R * N.pi * degrees * N.cos(N.pi*lat/180.)/180.

def m2deg(meters, lat=None):
    """Convert a distance in meters to degrees

    - **meters**: Distance in meters
    - *lat*: optional latitude, defaults to 0.

    Return: Distance in degrees

    .. seealso::

        :func:`deg2m`
    """

    if lat is None:
        lat = 0.
    else:
        lat = N.array(lat)

    return meters / (R * N.pi * N.cos(N.pi*lat/180.)/180.)

def basic_proj(lon, lat, inverse=False):
    """Convert a position from degrees to meters like a geographic projection

    :Params:

        - **lon**: Longitude in degrees or meters.
        - **lat**: Latitude in degrees or meters.
        - **inverse**, optional: Invert the projection from meters to degrees.

    """
    func = m2deg if inverse else deg2m
    return func(lon, lat), func(lat)



def ms2bf(ms):
    """Convert from m/s to Beauforts (wind)

    - **ms**: Wind speed in m/s.

    .. seealso::

        http://fr.wikipedia.org/wiki/%C3%89chelle_de_Beaufort
    """
    ms *= 3.6 # km/h to m/s
    return int(round((ms**2/9)**(1/3.)))

def dms2deg(d,m=0,s=0):
    """Convert from degrees/minutes/seconds to degrees

    - **d**: degrees
    - *m*: minutes [default: 0]
    - *s*: seconds [default: 0]

    .. seealso::

        :func:`deg2dms`
    """
    return (d < 0. and -1. or 1.) * (abs(d) + (m + s / 60.) /60.)

def deg2dms(deg):
    """Convert from degrees to degrees/minutes/seconds

    - **deg**: degrees

    .. seealso::

        :func:`dms2deg`
    """
    n, deg = deg < 0., abs(deg)
    d = int(deg)
    fm = (deg-d)*60.
    m = int(fm)
    s = (fm-m)*60.
    return n and -d or d, m, s

def kel2degc(tk):
    """Convert from degrees Kelvin to degrees Celsius

    - **tk**: Temperature in Kelvin.

    .. seealso::

        http://fr.wikipedia.org/wiki/Kelvin
    """
    #tk = tk - 273.15 # Kelvin to DegC
    return tk - 273.15

def degc2kel(dc):
    """Convert from degrees Celsius to degrees Kelvin

    - **dc**: Temperature in degrees Celsius.

    .. seealso::

        http://fr.wikipedia.org/wiki/Kelvin
    """
    #dc = dc + 273.15 # DegC to Kelvin
    return dc + 273.15

#def unorm(units):
#    """Try to normalize some units"""
#    units = units.strip()
#    if units.lower()=='m.p.h.': return 'mph'
#    m = re.match(r'\b(f(?:ee)?t|k?m|k(?:no)?t|kt)([\/\.\s]+|)(h(?:r)?|s(?:ec(?:ond)?)?|)(\b.*)', units.lower())
#    if m is not None:
#        return ''.join(m.groups())
#    return units

def tometric(units, value=1.,  munits=['m',  'm/s']):
    """Try to convert units to metric system using :mod:`~unidata.udunits.udunits`

    :Return: a float or ``None`` if conversion failed.
    """
    if isinstance(munits, basestring):
        munits = [munits]
    for mu in munits:
        try:
            return uuconvert(value, units, mu)
        except:
            pass

def rad2deg(r):
    return MV2.fmod(180 * r / MV2.pi, 360)

def deg2rad(d):
    return MV2.pi * MV2.fmod(d, 360) / 180

def vect2mod(u, v):
    return MV2.sqrt(MV2.power(u, 2) + MV2.power(v, 2))

def vect2dir(u, v):
    return rad2deg(MV2.arctan2(v, u))

def vect2moddir(u, v):
    return vect2mod(u, v), vect2dir(u, v)

def moddir2vectx(m, d):
    return MV2.multiply(m, MV2.cos(deg2rad(d)))

def moddir2vecty(m, d):
    return MV2.multiply(m, MV2.sin(deg2rad(d)))

def moddir2vectxy(m, d):
    return moddir2vectx(m, d), moddir2vecty(m, d)

# 1 kilooctet (ko) = 10**3 octets = 1 000 octets
sizeunits = {
    'K':10**3, 'M':10**6, 'G':10**9,
    'T':10**12, 'P':10**15, 'E':10**18,
    'Z':10**21, 'Y':10**24,
}

# 1 kibioctet (Kio) = 2**10 octets = 1 024 octets
sisizeunits = {
    'K':2**10, 'M':2**20, 'G':2**30,
    'T':2**40, 'P':2**50, 'E':2**60,
    'Z':2**70, 'Y':2**80,
}

def strfsize(size, fmt=None, si=None):
    """Format a size in bytes using the appropriate unit multiplicator (Ko, Mo, Kio, Mio)

    :Params:

        - **size**: the size in bytes
        - **fmt**: the format to use, will receive the size and the unit as format arguments
                   (None will automatically use "%.3f %s" or "%d %s")
        - **si**: whether to use International System units (10^3, ...) or not (2**10, ...)

    :Return: a string
    """
    if fmt is None:
        fmt = '%.3f %s' if float(size) % 1 else '%d %s'
    sortsizedict = lambda sd: reversed(sorted(sd.items(), lambda a, b: cmp(a[1],b[1])))
    if si is None: units,usfx = sortsizedict(sisizeunits),'o' # naive usage
    else: units,usfx = (sortsizedict(sisizeunits),'io') if si else (sortsizedict(sizeunits),'o')
    size = float(size)
    for unit,thresh in units:
        if size >= thresh:
            return fmt%(size/thresh, unit) + usfx
    return fmt%(size) + usfx

_strpsizerex = re.compile(r'(?P<number>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s*(?P<unit>%s)?(?P<usfx>io|o)?'%('|'.join(sizeunits.keys())), re.IGNORECASE)

def strpsize(size, si=True):
    """Parse a size in Ko, Mo, Kio, Mio, ...

    :Params:

        - **size**: the size string
        - **si**: whether to use International System units (10^3, ...) or not (2**10, ...)

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
            if s == 'io': return n * sisizeunits[u]
            elif si: return n * sisizeunits[u]
            else: return n * sizeunits[u]
        else: return n
    raise ValueError('Cannot parse size: %s'%(size))


