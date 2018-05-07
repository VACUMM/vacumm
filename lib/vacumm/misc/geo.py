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

from .misc import numod
from .constants import EARTH_RADIUS


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
    a = Nm.arctan2(Nm.cos(lat0)*Nm.sin(lat1) - Nm.sin(lat0)*Nm.cos(lat1)*Nm.cos(lon1-lon0),
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


def get_great_circle_points(lon0, lat0, lon1, lat1, npts, degrees=True, radius=None):
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
