# -*- coding: utf8 -*-
"""
Geospatial utilities
"""
# Copyright or © or Copr. Actimar/IFREMER (2018-2018)
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
import numpy as N
import MV2
import cdms2

from vacumm import VACUMMError
from .misc import numod
from .constants import EARTH_RADIUS
from .units import deg2m


def get_great_circle_angle(lon0, lat0, lon1, lat1, degrees=True, radius=None):
    """Get the great circle angular distance"""
    if degrees:
        lon0 *= N.pi/180
        lat0 *= N.pi/180
        lon1 *= N.pi/180
        lat1 *= N.pi/180
    Nm = numod(lon0, lat0, lon1, lat1)
    a = Nm.arccos(Nm.sin(lon0) * Nm.sin(lon1) +
                  Nm.cos(lon0) * Nm.cos(lon1) * Nm.cos(lat0 - lat1))
    if degrees:
        a *= 180 / N.pi
    return a


def haversine(lon0, lat0, lon1, lat1, degrees=True, radius=None):
    """Compute the haversine distance for a known radius"""
    if radius:
        radius = EARTH_RADIUS
    if degrees:
        lon0 *= N.pi/180
        lat0 *= N.pi/180
        lon1 *= N.pi/180
        lat1 *= N.pi/180
    Nm = numod(lon0, lat0, lon1, lat1)
    a = Nm.sin((lat0-lat1)/2)**2 + Nm.cos(lat0) * \
        Nm.cos(lat1) * Nm.sin((lon0-lon1)/2)**2
    return EARTH_RADIUS * 2 * Nm.arcsin(Nm.sqrt(a))


def get_bearing(lon0, lat0, lon1, lat1, degrees=True):
    """Compute the bearing angle (forward azimuth)"""
    if degrees:
        lon0 *= N.pi/180
        lat0 *= N.pi/180
        lon1 *= N.pi/180
        lat1 *= N.pi/180
    Nm = numod(lon0, lat0, lon1, lat1)
    a = Nm.arctan2(Nm.cos(lat0)*Nm.sin(lat1) -
                   Nm.sin(lat0)*Nm.cos(lat1)*Nm.cos(lon1-lon0),
                   Nm.sin(lon1-lon0)*Nm.cos(lat1))
    if degrees:
        a *= 180 / N.pi
    return a


def beardist2loc(lon0, lat0, bearing, dist, degrees=True, radius=None):
    if radius:
        radius = EARTH_RADIUS
    if degrees:
        lon0 *= N.pi/180
        lat0 *= N.pi/180
        bearing *= N.pi/180
    Nm = numod(lon0, lat0, bearing, dist)
    a = dist / EARTH_RADIUS
    lat1 = Nm.asin(Nm.sin(lat0)*Nm.cos(a) + Nm.cos(lat0)
                   * Nm.sin(a)*Nm.cos(bearing))
    lon1 = lon0 + Nm.atan2(Nm.cos(a)-Nm.sin(lat0)*Nm.sin(lat1),
                           Nm.sin(bearing)*Nm.sin(a)*Nm.cos(lat0))
    if degrees:
        lon1 *= 180 / N.pi
        lat1 *= 180 / N.pi
    return lon1, lat1


def get_great_circle_points(lon0, lat0, lon1, lat1, npts, degrees=True,
                            radius=None):
    """Get the great circle angle"""
    if radius:
        radius = EARTH_RADIUS
    if degrees:
        lon0 *= N.pi/180
        lat0 *= N.pi/180
        lon1 *= N.pi/180
        lat1 *= N.pi/180
    a = get_great_circle_angle(lon0, lat0, lon1, lat1, degrees=False)
    dists = N.reshape(N.linspace(0, 1, npts), N.shape(a)[::-1]+(1, )).T
    dists *= a * EARTH_RADIUS
    bearing = get_bearing(lon0, lat0, lon1, lat1, degrees=False)
    lon1, lat1 = beardist2loc(
        lon0, lat0, bearing, dists, degrees=False, radius=radius)
    if degrees:
        lon1 *= 180 / N.pi
        lat1 *= 180 / N.pi
    return lon1, lat1


def get_distances(xxa, yya, xxb=None, yyb=None, mode='haversine',
                  pairwise=False, geo=False):
    """Find the distances (in m) between a series of points

    Parameters
    ----------
    xxa:
        X coordinate of the first series
    yya:
        Y //
    xxb:
        X coordinate of the second series
    yyb:
        Y //
    mode: optional
        distance computation mode

        - ``None``: use ``"harversine"`` if longitude and latitude axes
          else ``"direct"``
        - ``"simple"`` or ``"euclidian"`` or ``"meters"``:
          simple euclidian distance with no coordinate
          tranformation
        - ``"harversine"`` or ``"sphere"`` or ``"degrees"``:
          great circle distance in meters from
          coordinates in degrees (see :func:`haversine`)
        - ``"deg2m"``: euclidian distance with coordinates converted
          from degrees to meters using :func:`~vacumm.misc.phys.units.deg2m`.
        - A callable object like a function that directly compute distance
          from coordinates::

                mode(xxa, yya, xxb, yyb)

    pairwise: optional
        The distances between A and B points are
        computed pairwise. There must have the same number of A and B points.

    Return
    ------
    array
        Distances as an ``(nb,na)`` array if pairwise if false,
        else a ``(na,)`` array.
    """
    # Mode
    if not callable(mode):
        if mode is None or geo:
            mode = "haversine"
        mode = str(mode).lower()
        valid_modes = ('simple', 'haversine', 'sphere', 'deg2m',
                       'degrees', 'euclidian', 'meters')
        if mode not in valid_modes:
            raise VACUMMError('Wrong mode ({}). '
                              'Please choose one of: {}'.format(
                                      mode, ', '.join(valid_modes)))
        if mode in ['sphere', "degrees"]:
            mode = 'harversine'
        elif mode in ["meters", "euclidian"]:
            mode = 'simple'

    # With what
    if xxb is None:
        xxb = xxa
    if yyb is None:
        yyb = yya

    # Numerical types
    Nma = numod(xxa, yya)
    Nmb = numod(xxb, xxb)
    nma = None if N.isscalar(xxa) and N.isscalar(yya) else Nma
    nmb = None if N.isscalar(xxb) and N.isscalar(yyb) else Nmb
    if MV2 is Nma or MV2 is Nmb:
        Nm = MV2
    elif N.ma is Nma or N.ma is Nmb:
        Nm = N.ma
    else:
        Nm = N

    # Make sur to have arrays
    oldshapea = N.shape(xxa) if N.ndim(xxa) >= N.ndim(yya) else N.shape(yya)
    oldshapeb = N.shape(xxb) if N.ndim(xxb) >= N.ndim(yyb) else N.shape(yyb)
    if cdms2.isVariable(xxa):
        xxa = xxa.asma()
        yya = yya.asma()
    else:
        xxa = N.atleast_1d(xxa)
        yya = N.atleast_1d(yya)
    if cdms2.isVariable(xxb):
        xxb = xxb.asma()
        yyb = yyb.asma()
    else:
        xxb = N.atleast_1d(xxb)
        yyb = N.atleast_1d(yyb)

    # Reshape them
    xxa = xxa.ravel().astype('d')
    yya = yya.ravel().astype('d')
    xxb = xxb.ravel().astype('d')
    yyb = yyb.ravel().astype('d')
    if not pairwise:
        xxa, xxb = N.meshgrid(xxa, xxb)
        yya, yyb = N.meshgrid(yya, yyb)

    # Compute it
    if callable(mode):
        dist = mode(xxa, yya, xxb, yyb)
    else:
        if mode == 'deg2m':
            xxa = deg2m(xxa, lat=yya)
            yya = deg2m(yya)
            xxb = deg2m(xxb, lat=yyb)
            yyb = deg2m(yyb)
        if mode == 'haversine':
            dist = haversine(xxa, yya, xxb, yyb, degrees=True)
        else:
            dist = Nm.sqrt((xxa-xxb)**2+(yya-yyb)**2)

    # Reform
    if nma or nmb:
        if pairwise:
            dist = dist.reshape(oldshapea or oldshapeb)
        else:
            dist = dist.reshape(oldshapeb + oldshapea)
    else:
        dist = float(dist)
    return dist
