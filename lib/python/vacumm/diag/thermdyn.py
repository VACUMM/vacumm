# -*- coding: utf8 -*-
"""Diagnostics about thermodynamics"""
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

from warnings import warn
import numpy as N, cdms2, MV2
from genutil.grower import grower
from vacumm import VACUMMError
from vacumm.data.cf import format_var, match_var
from vacumm.misc.grid import get_zdim, get_axis_slices, get_grid, set_grid
from vacumm.misc.axes import islat
from vacumm.misc import N_choose, grow_depth, grow_lat
from vacumm.misc.misc import grow_variables

try:
    from seawater import pres as sw_pres, dens as sw_dens, dens0 as sw_dens0, \
        dpth as sw_depth
except:
    try:
        from seawater.csiro import pres as sw_pres, dens as sw_dens, dens0 as sw_dens0, \
            depth as sw_depth
    except:
        warn('Failed to import seawater functions')


def density(temp, sal, depth=None, lat=None, potential=False,
    getdepth=False, getlat=False, format_axes=False):
    """Compute density from temperature, salinity and depth (and latitude)


    :Params:

        - **temp**: Insitu or potential temperature.
        - **sal**: Salinity.
        - **depth**, optional: Depth at temperature and salinty points.
          Assumed to be 0 if not found.
        - **lat**, optional: Latitude. Error when not found.
        - **potential**, optional: True to get the potential density (at atmospheric
          pressure).

    :Algo:

        >>> pressure = seawater.csiro.pres(depth, lat)
        >>> density = seawater.csiro.dens(sal, temp, depth)
    """

    # Compute
    if not potential and depth is not False: # In-situ

        # Get depth and latitude
        lat = grow_lat(temp, lat, mode='raise', getvar=False)
        if lat is None: raise VACUMMError('No latitude found for density')
        depth = grow_depth(temp, depth, mode='raise', getvar=False)
        if N.abs(depth.max())<N.abs(depth.min()): # positive
            depth = -depth
        if (depth.asma()<0).any():
            depth = depth-depth.min() # top=0

        # Get density
        pres = sw_pres(depth, lat)
        dens = sw_dens(sal, temp, pres) ; del pres

    else: # Potential

        dens = sw_dens0(sal, temp)
        getdepth = getlat = False

    # Format
    dens.setAxisList(temp.getAxisList())
    set_grid(dens, get_grid(temp))
    format_var(dens, 'dens', format_axes=format_axes)

    # Out
    if not getdepth and not getlat: return dens
    dens = dens,
    if getdepth: dens += depth,
    if getlat: dens += lat,
    return dens


def mixed_layer_depth(data, depth=None, lat=None, zaxis=None,
    mode=None, deltatemp=.2, deltadens=.03, kzmax=0.0005,
    potential=True, format_axes=False):
    """Get mixed layer depth from temperature and salinity

    :Params:

        - **temp**: Insitu or potential temperature.
        - **sal**: Salinity.
        - **depth**, optional: Depth at temperature and salinty points.
        - **lat**, optional: Latitude.
        - **mode**, optional: ``"deltatemp"``, ``"deltadens"``, ``"kz"``
          or ``"twolayers"``


    :Raise: :class:`~vacumm.VACUMMError` if can't get depth (and latitude for density).
    """

    # TODO: positive up

    # Inspection
    if isinstance(data, tuple): # data = temp,sal

        temp, sal=data

        # Get density
        if mode!='deltatemp':

            res = density(temp, sal, depth=depth, lat=lat,
                format_axes=False, potential=potential, getdepth=True)
            if isinstance(res, tuple):
                dens, depth = res
            else:
                dens = res
            dens = dens.asma()
            if mode is None:
                mode = 'deltadens'

        else:

            temp = data[0]

        # Check mode
        if mode == 'kz':
            warn("Switching MLD computation mode to 'deltadens'")
            mode = "deltadens"

    elif match_var(data, 'temp', mode='nslu'):

        if mode is not None and mode!='deltatemp':
            warn("Switching MLD computation mode to 'deltatemp'")
        mode = 'deltatemp'
        temp = data

    elif match_var(data, 'dens', mode='nslu'):

        if mode in ['kz', 'deltatemp']:
            warn("Switching MLD computation mode to 'deltadens'")
            mode = None
        if mode is None:
            mode = "deltadens"
        dens = data

    elif match_var(data, 'kz', mode='nslu'):

        if mode is None:
            mode = "kz"
        if mode != "kz":
            warn("Switching MLD computation mode to 'kz'")
        kz = data

    else:

        if mode in ['deltadens', 'twolayers']:
            dens = data
        elif mode == "deltatemp":
            temp = data
        elif mode == "kz":
            kz = data
        elif mode is not None:
            raise VACUMMError("Invalid MLD computation mode : '%s'"%mode)
        else:
            raise VACUMMError("Can't guess MLD computation mode")

        temp = delta

    # Find Z dim
    data0 = data[0] if isinstance(data, tuple) else data
    depth = grow_depth(data0, depth, mode='raise', getvar=False)
    zaxis = get_zdim(data0, axis=zaxis)
    if zaxis is None:
        raise VACUMMError("Can't guess zaxis")
    slices = get_axis_slices(data0, zaxis)

    # Init MLD
    axes = data0.getAxisList()
    del axes[zaxis]
    mld = MV2.array(data0.asma()[slices['first']], copy=1, axes=axes, copyaxes=False)
    set_grid(mld, get_grid(data0))
    format_var(mld, 'mld', format_axes=format_axes)
    mld[:] = MV2.masked

    # Two-layers
    if mode=='twolayers':

        densbot = dens[slices['first']]
        denstop = dens[slices['last']]
        del dens
        H = 1.5*depth[slices['first']] - 0.5*depth[slices['firstp1']]
        H = -1.5*depth[slices['last']] + 0.5*depth[slices['lastm1']]
        mld[:] = -H*(densbot-denstop)/(densbot-denstop)
        del H

    elif mode=='deltadens':

        denscrit = dens[slices['last']]+deltadens
        mld[:] = -_val2z_(dens, depth, denscrit, zaxis, -1)
        del dens

    elif mode=='deltatemp':

        tempcrit = temp[slices['last']]-deltatemp
        mld[:] = -_val2z_(temp, depth, tempcrit, zaxis, 1)

    elif mode=='kz':

        mld[:] = -_valmin2z_(kz, depth, kzmax, zaxis, 1)

    else:

        raise VACUMMError("Invalid mode for computing MLD (%s)."%mode +
            "Please choose one of: deltadens, twolayers")

    # Mask zeros
    mld[:] = MV2.masked_values(mld, 0., copy=0)

    return mld



def _val2z_(var, dep, varref, axis, monosign):
    """Compute depth for a given reference data value"""

    # Inits with Z axis first
    var = var.asma() if cdms2.isVariable(var) else N.ma.asarray(var)
    var = var.copy()
    dep = dep.asma() if cdms2.isVariable(dep) else N.ma.asarray(dep)
    dep = N.ma.masked_where(var.mask, dep, copy=False)
    varref = varref.asma() if cdms2.isVariable(varref) else N.ma.asarray(varref)
    if var.ndim-1==0:
        var = var.reshape(-1,1)
        dep = dep.reshape(-1,1)
        varref = varref.reshape(1)
        axis = 0
    rvar = N.rollaxis(var, axis) if axis else var
    rdep = N.rollaxis(dep, axis) if axis else dep
    rdepref = rdep[0].copy()
    rvar = rvar.copy()
    nz = rvar.shape[0]   # nb of elements along the first axis (z in general)
    maskland = N.ma.getmaskarray(rvar)[-1]

    # Monotonic
    rvar *= monosign
    varref *= monosign
    dvar = N.ma.diff(rvar, axis=0)

    # Find z index
    dvarref = rvar-varref
    iiz0 = N.arange(nz-1, dtype='i')
    iiz0 = N.ma.resize(iiz0, rvar.shape[rvar.ndim-1:0:-1]+(nz-1, )).T
    iiz0 = N.ma.masked_where(N.ma.diff(N.sign(dvarref), axis=0)<2, iiz0, copy=0)
    iz0 = iiz0.max(axis=0)
    del iiz0
    maskdz = iz0.mask
    iz0 = iz0.filled(0.)


    # Interpolation
    dvar0 = N_choose(iz0, dvar)
    dvar0[maskdz] = 1.
    w1 = -N_choose(iz0, dvarref)/dvar0
    w1.set_fill_value(0.)
    rdep.set_fill_value(0.)
    rdepref = N_choose(iz0+1, rdep)*w1
    rdepref += N_choose(iz0, rdep)*(1-w1)

    # Masking
    rdepref[maskland] = N.ma.masked
    rdepref = N.ma.where(~maskland&maskdz, rdep.min(axis=0), rdepref)

    return rdepref


def _valmin2z_(var, dep, varmax, axis, monosign):
    """Compute depth at which data value at lower than a reference value"""

    # Inits with Z axis first
    var = var.asma() if cdms2.isVariable(var) else N.ma.asarray(var)
    var = var.copy()
    dep = dep.asma() if cdms2.isVariable(dep) else N.ma.asarray(dep)
    if var.ndim-1==0:
        var = var.reshape(-1,1)
        dep = dep.reshape(-1,1)
        axis = 0
    rvar = N.rollaxis(var, axis) if axis else var
    rdep = N.rollaxis(dep, axis) if axis else dep
    rvar = rvar.copy()
    nz = rvar.shape[0]
    maskland = N.ma.getmaskarray(rvar)[-1]

    # Monotonic (greater at surface)
    rvar *= monosign

    # Find z index
    iiz = N.arange(nz, dtype='i')
    iiz = N.ma.resize(iiz, rvar.shape[rvar.ndim-1:0:-1]+(nz, )).T
    iiz = N.ma.masked_where(rvar>varmax, iiz, copy=0)
    iz0 = iiz.max(axis=0)
    del iiz
    maskz = iz0.mask
    iz0 = iz0.filled(0.)

    # Values
    rdepref = N_choose(iz0, rdep)

    # Masking
    rdepref[maskland] = N.ma.masked
    rdepref = N.ma.where(~maskland&maskz, rdep.min(axis=0), rdepref)

    return rdepref

if __name__=='__main__':



#    print 'hoho'
    nx = 4
    ny = 2
    nz = 6
    depth = N.ma.arange(nz*nx*ny**1.).reshape(nz, ny, nx)
    depth = depth-depth.max()
    # partie ci-dessus inutile en fait

#    var = N.round(N.sqrt(-depth)*10)
#    var[:,1, 2] = N.ma.masked
#
#    varref = var[-1]+44#15.7
#    varref[:] = 23
#    varref[:] = 67.5
#
#    var = var[::-1]
#    varref = var[-1]*0+67.5
#
    # tableaux 1D (29 elements)
    var = N.ma.array([13.1168556213, 13.120059967, 13.1211280823, 13.1211280823, 13.1211280823
, 13.1211280823, 13.1179237366, 13.1136512756, 13.1093788147, 13.1083106995
, 13.1072416306, 13.1019010544, 13.066652298, 13.0527667999, 13.0784015656
, 13.1136512756, 13.1905574799, 13.2461013794, 13.2503738403, 13.2236700058
, 13.2236700058, 13.5729541779, 14.4798116684, 15.1249732971, 15.0694293976
, 15.0683612823, 15.0672931671, 15.0651569366, 15.0608844757])
    depth = N.array([-2340.67651367, -2240.55786133, -2134.18212891, -2021.54858398
, -1908.9152832, -1796.28186035, -1683.64868164, -1571.01525879, -1458.38183594
, -1345.74865723, -1233.11523438, -1120.48181152, -1007.84857178, -895.215148926
, -782.581848145, -669.948547363, -563.572570801, -463.454071045, -369.592926025
, -284.4921875, -209.403320312, -150.583679199, -106.781829834, -75.4947891235
, -50.4651565552, -27.9384918213, -12.9207115173, -5.41182279587, -1.65737783909])

    varref = var[-1]-0.2
    print varref, var[-1]

    var = N.ma.array([13.1157875061, 13.1157875061, 13.1157875061, 13.1157875061, 13.1157875061
, 13.1168556213, 13.1168556213, 13.1168556213, 13.1168556213, 13.1179237366
, 13.1189918518, 13.1275367737, 13.1414222717, 13.1168556213, 13.1019010544
, 13.1798763275, 13.2610549927, 13.2589187622, 13.2450332642, 13.2428970337
, 13.2450332642, 13.5654773712, 14.426404953, 15.0363168716, 14.9070711136
, 14.9401836395, 14.9700918198, 14.9668874741, 14.9604787827])
    depth = N.ma.array([-2478.62475586, -2372.61083984, -2259.97119141, -2140.70556641
, -2021.44006348, -1902.17443848, -1782.90893555, -1663.64331055, -1544.37768555
, -1425.11218262, -1305.84655762, -1186.58093262, -1067.31530762, -948.049682617
, -828.784118652, -709.518554688, -596.87878418, -490.864959717, -391.476959229
, -301.365203857, -221.854782104, -159.571640015, -113.190582275, -80.0612487793
, -53.5577850342, -29.7046699524, -13.8025913239, -5.85155200958, -1.87603235245])
    varref = var[-1]-0.2
    print varref, var[-1]

    depref = _val2z_(var, depth, varref, 0, 1)

#    print 'var',var[:,0]
#    print 'varref', varref[0]
#    print 'dep',depth[:,0]
#    print 'res',depref[0]
    print depref
