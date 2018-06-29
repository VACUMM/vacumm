# -*- coding: utf8 -*-
"""Statistical tools"""

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

import os

import numpy as N,numpy.ma as MA,MV2
import cdms2
import cdtime
from genutil.statistics import percentiles

from vacumm import VACUMMError
from vacumm.misc.axes import isaxis, create_time
from vacumm.misc.grid import get_grid, set_grid
from vacumm.misc import MV2_axisConcatenate, set_atts, cp_atts
from vacumm.misc.io import netcdf4
from vacumm.misc.atime import comptime, reltime

__all__ = ['corr_proba', 'ensrank', 'qtmin', 'qtmax', 'qtminmax',
    'StatAccum', 'StatAccumError']
__all__.sort()

def corr_proba(r, ndata, ndataset=2, dof=False):
    """Probability of rejecting correlations

    - **r**: Correlation coefficient
    - **ndata**: Number of records use for correlations
    - **ndataset**, optional:  Number of datasets (1 for autocorrelations, else 2) [default: 2]

    .. todo::

        This must be rewritten using :mod:`scipy.stats`
    """
    # TODO: use scipy for betai and _gamma?
    from genutil.salstat import betai,_gammaln

    # Basic tests
    ndata = MA.masked_equal(ndata,0,copy=0)
    r = MV2.masked_where(MA.equal(MA.absolute(r),1.),r,copy=0)

    # Degree of freedom
    if dof:
        df = ndata
    else:
        df = ndata-2-ndataset

    # Advanced test: prevent extreme values by locally decreasing the dof
    reduc = N.ones(r.shape)
    z = None
    while z is None or MA.count(MA.masked_greater(z,-600.)):
        if z is not None:
            imax = MA.argmin(z.ravel())
            reduc.flat[imax] += 1
        dfr = df/reduc
        t = r*MV2.sqrt(dfr/((1.0-r)* (1.0+r)))
        a = 0.5*dfr
        b = 0.5
        x = df/(dfr+t**2)
        z = _gammaln(a+b)-_gammaln(a)-_gammaln(b)+a*MA.log(x)+b*MA.log(1.0-x)

    # Perfom the test and format the variable
    prob = MV2.masked_array(betai(a,b,x),axes=r.getAxisList())*100
    prob.id = 'corr_proba' ; prob.name = prob.id
    prob.long_name = 'Probability of rejection'
    prob.units = '%'

    return prob


def ensrank(obs, ens, gethist=False, getnrz=False, centered=False):
    """Compute the rank of a reference (observations) in an ensemble

    If ``nens`` is the size of the ensemble,
    the rank may go from ``0`` to ``nens``.

    - **obs**, (...): Observation array
    - **ens**, (nens,...): Ensemble array
    - **getrank**, optional: If True, return also the rank histogram
    - **getnrz**, optional: If True, return also non recovery zones
    - **centered**, optional: Center the ensemble data on the observations
      before computing the rank
    """
    # Check
    ens = MV2.asarray(ens)
    obs = MV2.asarray(obs)
    nens = ens.shape[0]
    assert nens>3, 'Your ensemble is too small (nens=%i)'%nens

    # Inits
    rank = obs.clone()
    rank.id = 'rank'
    rank.long_name = 'Rank'
    if hasattr(rank, 'units'): del rank.units
    nrz = obs.clone()
    nrz.id = 'nrz'
    nrz.long_name = 'Non recovery zones'
    nrz[:] = MV2.masked

    # Sort ensemble at each pixel and compute the min biais
    bias = MV2.sort(ens, axis=0)-obs
    if centered:
        bias -= bias.mean(axis=0)
    nrz[:] = MV2.where((bias[-1]<0.).astype('b'), bias[-1], nrz)
    nrz[:] = MV2.where((bias[0]>0.).astype('b'), bias[0], nrz)


    # Index of member  = rank
    rank[:] = (bias<0).astype('i').sum(axis=0)
    del bias

    ## Apply observation mask
    #nrz[:] = MV2.masked_where(MV2.getmask(obs), nrz, copy=0)
    #rank[:] = MV2.masked_where(MV2.getmask(obs), rank, copy=0)

    # Rank histogram
    if gethist:
        if N.ma.isMA(rank): rank = rank.compressed()
        hist = N.histogram(rank, N.arange(0,nens+1))

    # Out
    ret = rank,
    if gethist:
        ret += hist,
    if getnrz:
        ret += nrz,
    return ret[0] if len(ret)==1 else ret


def qtminmax(var, qt=5.):
    """Get min and max using quantiles

    This is useful for very asymmetic distributions

    :Params:

        - **var**: A numeric array.
        - **qt**, optional: Percentile used to define min and max

            - Single number: ``qt`` -> ``[qt, 100-qt]``
            - Two elements: used directly.
    """
    if N.ma.isMA(var): var = var.compressed()
    if not hasattr(qt, '__len__'):
        qt = [qt, 100-qt]
    qt = N.sort(N.clip(N.asarray(qt), 0, 100.))
    vmin, vmax = percentiles(var.ravel(), qt)
    return vmin, vmax

def qtmin(var, qt=5.):
    """Get min max using quantiles

    This is useful for very asymmetic distributions

    :Params:

        - **var**: A numeric array.
        - **qt**, optional: Percentile used to define min
    """
    if N.ma.isMA(var): var = var.compressed()
    qt = N.clip([qt], 0, 100.)
    vmin = percentiles(var.ravel(), qt)
    return vmin

def qtmax(var, qt=95.):
    """Get max using quantiles

    This is useful for very asymmetic distributions

    :Params:

        - **var**: A numeric array.
        - **qt**, optional: Percentile used to define max
    """
    if N.ma.isMA(var): var = var.compressed()
    qt = N.clip([qt], 0, 100.)
    vmax = percentiles(var.ravel(), qt)
    return vmax


class StatAccumError(VACUMMError):
    pass

class StatAccum(object):
    """Statistics accumulator

    :Generic params:

        - **t/sall**: Perform all the statistics by default.
        - **withtime**, optional: Input does not contain a time dimension.

    :Single:
        - **tcount/scount**: Number of observations taken into account
        - **tavail/savail**: Percentage of available observations
        - **tmean/smean**: Temporal (t) / Spatial (s) average
        - **tstd/sstd**: Temporal (t) / Spatial (s) std

    :Dual:
        - **tall**: Perform all the following statistics by default.
        - **tbias/sbias**: Temporal (t) / Spatial (s) bias
        - **trms/srms**: Temporal (t) / Spatial (s) RMS
        - **tcrms/scrms**: Temporal (t) / Spatial (s) centered RMS
        - **tcorr/scorr**: Temporal (t) / Spatial (s) Correlation

    :Example:


        >>> sa = StatAccum(tbias=True, srms=True)
        >>> sa += ssta1, sstb1
        >>> sa += ssta2, sstb2
        >>> tbias = sa.get_tbias()
        >>> srms = sa.get_srms()

    """

    single_stats = 'mean', 'std', 'hist', 'min', 'max'
    dual_stats = 'bias', 'rms', 'crms', 'corr', 'avail', 'count'
    all_stats = single_stats + dual_stats
    _single_accums = 'sum', 'sqr', 'min', 'max', 'hist'
    _dual_accums = 'prod',
    default_restart_file = 'stataccum_restart.nc'
    def __init__(self,
        tall=False, tavail=None, tmean=None, tstd=None, smax=None, smin=None,
        tbias=None, trms=None, tcrms=None, tcorr=None, shist=None, tcount=None,
        sall=False, savail=None, smean=None, sstd=None, tmax=None, tmin=None,
        sbias=None, srms=None, scrms=None, scorr=None, thist=None, scount=None,
        bins=None, restart_file=None, restart=False,
        withtime=None):

        # Restart?
        if restart:
            if restart_file is None and isinstance(restart, basestring):
                restart_file = restart
            self.restart_file = restart_file
            self.load()
        else:

            for st in 'st':
                activate = False
                for stype in self.all_stats:
                    value = eval(st+stype)
                    if value is None: value = eval(st+'all')
                    setattr(self, st+stype, value)
                    activate |= value
                setattr(self, st+'stats', activate)
            self.iterindex = 0
            self.lasttime = None
            self.withtime = withtime
            if self.thist or self.shist:
                if bins is None:
                    if (self.thist and not self.tall) or (self.shist and not self.sall):
                        raise StatAccumError(
                            "You must set the 'bins' keyword to make histograms")
                    if self.thist and self.tall:
                        self.thist = False
                    if self.shist and self.sall:
                        self.shist = False
                self._baxis = cdms2.createAxis(0.5*(bins[:-1]+bins[1:]), id='hbin')
                self._baxis.long_name = 'Histogram bin centers'
                self._baxis.setBounds(bins)
            else:
                self._baxis = None
            self.bins = bins
            self.nitems = 0
            self.dual = None
            self.restart_file = None
            self.ns = self.nt = -1

    def __len__(self):
        return self.iterindex

    def append(self, *items):
        """Append data to the accumulator"""
        # Initialization
        if self.iterindex==0:

            # Two variables?
            self.nitems = len(items)
            self.dual = self.nitems>1

            # Do we have time in input variables
            if self.withtime is None:
                self.withtime = items[0].getTime() is not None

            # Flags
            if not self.dual:
                self.tbias = self.trms = self.tcrms = self.tcorr = False
                self.sbias = self.srms = self.scrms = self.scorr = False
            self.tstats = (self.tmean or self.tstd or self.tbias or self.trms or
                self.tcrms or self.tavail or self.tmin or self.tmax or self.thist or
                self.tcount)
            self.sstats = (self.smean or self.sstd or self.sbias or self.srms or
                self.scrms or self.savail or self.smin or self.smax or self.shist or
                self.scount)
            self.tsum = self.tmean or self.tbias or self.tstd or self.tcrms or self.tcorr
            self.tsqr = self.tstd or self.trms or self.tcrms or self.tcorr
            self.tprod = self.trms or self.tcrms or self.tcorr
            self.ssum = self.smean or self.sbias or self.sstd or self.scrms or self.scorr
            self.ssqr = self.sstd or self.srms or self.scrms or self.scorr
            self.sprod = self.srms or self.scrms or self.scorr
            if self.shist or self.thist:
                self.nbins = self.bins.shape[0]-1

            # Check spatial statistics
            if self.sstats and items[0].ndim<2:
                self.sstats = self.smean = self.sstd = self.sbias = self.srms = \
                    self.scrms = self.scorr = self.savail = self.smin = self.smax = \
                    self.shist = self.scount = False

            # Init template for output temporal statistics
            if self.tstats:
                self._ttemplates = ()
                self._thtemplates = None
                if self.thist:
                    self._thtemplates = ()
                selspace = 0 if self.withtime else slice(None)
                for item in items:
                    tpl = item[selspace].clone()
                    for att in tpl.attributes:
                        delattr(tpl, att)
                    self._ttemplates += tpl,
                    if self.thist:
                        htpl = self._template_t2ht_(tpl)
                        self._thtemplates += htpl,

                self._tbase = items[0].filled()[selspace]*0.
                self._tbase.shape = -1
                if self.thist:
                    self._thbase = N.resize(self._tbase.astype('l'),
                        (self.nbins, self._tbase.size))

            # Save some attributes
            self._atts = ()
            for var in items:
                attr = {}
                for att in 'id', 'units', 'long_name', 'standard_name':
                    if hasattr(var, att):
                        attr[att] = getattr(var, att)
                self._atts += attr,


            # Init accumulated variables

            # - temporal statistics
            if self.tstats:
                self._nts = []
                if self.tsum:
                    self._tsum = self._tbase.copy(),
                    if self.dual: self._tsum += self._tbase.copy(),
                # - square
                if self.tsqr:
                    self._tsqr = self._tbase.copy(),
                    if self.dual: self._tsqr += self._tbase.copy(),
                # - product
                if self.tprod: self._tprod = self._tbase.copy()
                # - count
                self._tcount = self._tbase.copy().astype('l')
                self.nt = 0
                # - min/max
                if self.tmin:
                    self._tmin = self._tbase.copy()+N.inf,
                    if self.dual: self._tmin += self._tbase.copy()+N.inf,
                if self.tmax:
                    self._tmax = self._tbase.copy()-N.inf,
                    if self.dual: self._tmax += self._tbase.copy()-N.inf,
                # - histograms
                if self.thist:
                    self._thist = self._thbase.copy(),
                    if self.dual:
                        self._thist += self._thbase.copy(),
            # - spatial statistics
            if self.sstats:
                self._sstats = {}
#                for name in self.all_stats:
#                    if not getattr(self, 's'+name): continue
#                    self._sstats[name] = []
#                    setattr(self, '_s'+name, [])
                self._stimes = tuple([[] for i in xrange(len(items))]) if self.withtime else None
                self._scount = []
            self.ns = items[0].size/items[0].shape[0]

        # Checks
        if not self.tstats and not self.sstats:
            raise StatAccumError("Cant' accumulate statistics since they are all are OFF.")
        if self.dual:
            if len(items)!=2:
                raise StatAccumError("You must always provide two variables for accumulating stats")
            if items[0].shape!=items[1].shape:
                raise StatAccumError("Compared variables must have the same shape: %s != %s"%(
                    items[0].shape, items[1].shape))

        # As 2D (time-space) zero-filled arrays
        nt = items[0].shape[0] if self.withtime else 1
        self._nts.append(nt)
        item0 = items[0].asma().reshape(nt, -1)
        mask0 = N.ma.getmaskarray(item0)
        if self.dual:
            item1 = items[1].asma().reshape(nt, -1)
            mask1 = N.ma.getmaskarray(item1)
            mask = mask0 | mask1
        else:
            mask = mask0
        item0[mask] = N.ma.masked
        if self.tmax or self.tmin or self.smax or self.smin:
            item0m = item0
        item0 = item0.filled(0.)
        if self.dual:
            item1[mask] = N.ma.masked
            if self.tmax or self.tmin or self.smax or self.smin:
                item1m = item1
            item1 = item1.filled(0.)

        # Accumulate

        # - temporal statistics
        if self.tstats:
            self._tcount += (~mask).astype('l').sum(axis=0)
            self.nt += items[0].shape[0]
            if self.tsum:
                self._tsum[0][:] += item0.sum(axis=0)
                if self.dual: self._tsum[1][:] += item1.sum(axis=0)
            if self.tsqr:
                self._tsqr[0][:] += (item0**2).sum(axis=0)
                if self.dual: self._tsqr[1][:] += (item1**2).sum(axis=0)
            if self.tprod:
                self._tprod[:] += (item0*item1).sum(axis=0)
            if self.tmin:
                vmin = item0m.min(axis=0)
                self._tmin[0][:] = N.ma.where(vmin<self._tmin[0], vmin, self._tmin[0])
                if self.dual:
                    vmin = item1m.min(axis=0)
                    self._tmin[1][:] = N.ma.where(vmin<self._tmin[1], vmin, self._tmin[1])
                del vmin
            if self.tmax:
                vmax = item0m.max(axis=0)
                self._tmax[0][:] = N.ma.where(vmax>self._tmax[0], vmax, self._tmax[0])
                if self.dual:
                    vmax = item1m.max(axis=0)
                    self._tmax[1][:] = N.ma.where(vmax>self._tmax[1], vmax, self._tmax[1])
                del vmax
            if self.thist:
                indices0 = N.digitize(item0.flat, self.bins)
                indices0[mask.ravel()] = 0
                if self.dual:
                    indices1 = N.digitize(item1.flat, self.bins)
                    indices1[mask.ravel()] = 0
                for ib in xrange(self.nbins):
                    valid0 = indices0.reshape(item0.shape)==ib+1
                    valid0 &= ~mask
                    self._thist[0][ib] += valid0.sum(axis=0)
                    if self.dual:
                        valid1 = indices1.reshape(item1.shape)==ib+1
                        valid1 &= ~mask
                        self._thist[1][ib] += valid1.sum(axis=0)
                del indices0, valid0
                if self.dual: del indices1, valid1


        # - spatial statistics
        if self.sstats:

            # Count
            scount = (~mask).astype('l').sum(axis=1)
#            if self.savail:
            self._scount.append(scount)

            # Base
            if self.ssum:
                ssum = item0.sum(axis=1),
                if self.dual: ssum += item1.sum(axis=1),
            if self.ssqr:
                ssqr = (item0**2).sum(axis=1),
                if self.dual:
                    ssqr += (item1**2).sum(axis=1),
            if self.sprod:
                sprod = (item0*item1).sum(axis=1)
            if self.smin:
                smin = item0m.min(axis=1),
                if self.dual:
                    smin += item1m.min(axis=1),
            if self.smax:
                smax = item0m.max(axis=1),
                if self.dual:
                    smax += item1m.max(axis=1),
            if self.shist:
                indices0 = N.digitize(item0.flat, self.bins)
                indices0[mask.ravel()] = 0
                if self.dual:
                    indices1 = N.digitize(item1.flat, self.bins)
                    indices1[mask.ravel()] = 0

                shist = N.zeros(item0.shape[:1]+(self.nbins, ), 'l'),
                if self.dual: shist += shist[0].copy(),

                for ib in xrange(self.nbins):
                    valid0 = indices0.reshape(item0.shape)==ib+1
                    shist[0][:, ib] += valid0.sum(axis=1)
                    if self.dual:
                        valid1 = indices1.reshape(item1.shape)==ib+1
                        shist[1][:, ib] += valid1.sum(axis=1)


            # Stats
            sstats = self._base2stats_(self.ns, scount, ssum if self.ssum else None,
                ssqr if self.ssqr else None, sprod if self.sprod else None,
                smin if self.smin else None, smax if self.smax else None,
                shist if self.shist else None,
                self.smean, self.sstd, self.sbias, self.srms, self.scrms,
                self.scorr, self.savail, self.scount,
                self.smin, self.smax, self.shist)
            for key, val in sstats.iteritems():
#                self._sstats[key].append(val)
#                getattr(self, '_s'+key).append(val)
                if isinstance(val, tuple):
                    self._sstats.setdefault(key, ([], []))
                    for i, v in enumerate(val):
                        self._sstats[key][i].append(v)
                else:
                    self._sstats.setdefault(key, [])
                    self._sstats[key].append(val)

            # Time axes (for concatenation)
            if self.withtime:
                for i, item in enumerate(items):
                    self._stimes[i].append(item.getAxis(0))

        # Save iterative status
        self.iterindex += 1
        if self.withtime:
            self.lasttime = items[0].getTime().asComponentTime()[-1]


    def __iadd__(self, items):
        if isinstance(items, tuple):
            self.append(*items)
        else:
            self.append(items)
        return self

    def isempty(self):
        return self.iterindex==0

    def _template_t2ht_(self, tpl):
        htpl = MV2.resize(tpl.astype('l'), (self.nbins,)+tpl.shape)
        htpl.setAxisList([self._baxis]+tpl.getAxisList())
        set_grid(htpl, get_grid(tpl))
        cp_atts(tpl, htpl)
        return htpl

    def _asarray_(self, a):
        """Merge a list of arrays along time"""
        if isinstance(a, list):
            return N.concatenate(a, axis=0)
        return a

    def _aslist_(self, a):
        """Split an array along time"""
        return N.split(a, N.cumsum(self._nts[:-1]), axis=0)

    def _dump_array_(self, f, var, id, axis):
        """Dump an array or a list of them in a nctdf file"""
        istime = isinstance(var, list)
        ishist = 'hist' in id
        var = self._asarray_(var)
        var = MV2.array(var, copy=0, id=id)
        var.setAxis(int(not istime and ishist), axis)
        if ishist:
            var.setAxis(istime, self._baxis)
        if var.dtype.char in 'fd':
            var.set_fill_value(1e20)
        else:
            var.set_fill_value(-1)
        return f.write(var, extend=0)

    def _load_array_(self, f, id):
        var = f(id)
        istime = 't' in var.getOrder()
        var = var.filled()
        if istime:
            var = self._aslist_(var)
        return var

    def dump(self, restart_file=None):
        """Dump the current instance to a netcdf file"""

        # File
        netcdf4()
        if restart_file is None:
            restart_file = self.restart_file
        if restart_file is None:
            restart_file = self.default_restart_file
        self.restart_file = restart_file
        if os.path.exists(restart_file):
            os.remove(restart_file)
        f = cdms2.open(restart_file, 'w')

        # Config
        # - what was initially asked and some more
        for sname in self.all_stats + ('sum', 'sqr', 'prod', 'stats'):
            for st in 'st':
                value = getattr(self, st+sname)
                if value is None: continue
                if isinstance(value, bool): value = int(value)
                setattr(f, st+sname, value)
        # - current status
        f.iterindex = self.iterindex
        f.nitems = int(self.nitems)
        f.withtime = -1 if self.withtime is None else int(self.withtime)
        if self.withtime and self.lasttime:
            f.lasttime = str(self.lasttime)
        if self.bins is None:
            f.bin_edges = 0
        else:
            f.bin_edges = self.bins
        if f.nitems==0: # Still no data
            f.close()
            return
        # - already had some data
        f.dual = int(self.dual)
        f.ns = self.ns
        f.nt = self.nt
        f.nts = self._nts
        f.tstats = int(self.tstats)
        f.sstats = int(self.sstats)

        # Spatial statistics
        if self.sstats:

            # Time axes
            if self.withtime:
                taxes = ()
                for i, tt in enumerate(self._stimes):
                    taxis = MV2_axisConcatenate(tt)
                    taxis.stataccum_oldid = taxis.id
                    taxis.id = 't'+str(i)
                    taxes += taxis,
            else:
                taxes = [cdms2.createAxis(N.arange(self.nt), id='t')]

            # Count
            self._dump_array_(f, var=self._scount, id='scount', axis=taxes[0])

            # Other stats
            for key, stats in self._sstats.items():
                if isinstance(stats, tuple):
                    for i, stat in enumerate(stats):
                        self._dump_array_(f, var=stat, id='s%s%s'%(key, str(i)),
                            axis=taxes[i])
                else:
                    self._dump_array_(f, var=stats, id='s'+key, axis=taxes[0])

        # Temporal statistics
        if self.tstats:

            # Main axis
            axis = cdms2.createAxis(N.arange(self.ns), id='s')

            # Count
            self._dump_array_(f, var=self._tcount, id='tcount', axis=axis)

            # Other stats
            for key in self._dual_accums+self._single_accums:
                if not getattr(self, 't'+key): continue
                value = getattr(self, '_t'+key)
                if key in self._dual_accums:
                    self._dump_array_(f, var=value, id='t'+key, axis=axis)
                else:
                    for i, stat in enumerate(value):
                        self._dump_array_(f, var=stat, id='t%s%s'%(key, str(i)), axis=axis)


            # Templates
            ttemplates = [self._ttemplates]
#            if self.thist:
#                ttemplates.append(self._thtemplates)
            for ttpls in ttemplates:
                for i, ttpl in enumerate(ttpls):
#                    ttpl = ttpl.clone(copyData=0)
#                    for ia, axis in enumerate(ttpl.getAxisList()): #FIXME:grid cloning
#                        ttpl.setAxis(ia, axis.clone(copyData=0))
                    _add_id_prefix_(ttpl, 'var%i_'%i, exc=self._baxis)
                    f.write(ttpl)
                    _rm_id_prefix_(ttpl, 'var%i_'%i, exc=self._baxis)

        # Attributes
        for ivar, atts in enumerate(self._atts):
            var = MV2.zeros(1)
            set_atts(var, atts)
            var.id = 'var%i_atts'%ivar
            var.stataccum_id = atts['id']
            var.getAxis(0).id = 'var_att_axis'
            f.write(var)

#        self._ttemplates, self._thtemplates, self._atts, self._tbase

        f.close()
        return f.id


    def load(self, restart_file=None, iterindex=None, nowtime=None):
        """Load the current instance from a netcdf file

        :Params:

            - **restart_file**, optional: Netcdf restart file.
            - **iterindex**, optional: If given, the restart file is not loaded if
              ``iterindex`` is greater or equal to the file's ``iterindex`` attribute.
            - **nowtime**, optional: If given, the restart file is not loaded if
              ``nowtime`` is greater or equal to the file's ``lasttime`` attribute.
        """

        # File
        if restart_file is None:
            restart_file = self.restart_file
        if restart_file is None:
            restart_file = self.default_restart_file
        self.restart_file = restart_file
        f = cdms2.open(restart_file)

        # Config
        # - check status
        if iterindex is not None:
            self.iterindex = iterindex
        if hasattr(self, 'iterindex') and f.iterindex<=self.iterindex:
            return -1
        if nowtime is not None:
            self.lasttime = comptime(nowtime)
        if (hasattr(self, 'lasttime') and f.withtime>0 and self.lasttime
                and reltime(f.lasttime, 'hours since 2000').value <=
                reltime(self.lasttime, 'hours since 2000').value):
            return -1
        # - what was initially asked and some more
        for sname in self.all_stats + ('sum', 'sqr', 'prod', 'stats'):
            for st in 'st':
                if not hasattr(f, st+sname): continue
                value = getattr(f, st+sname)
                setattr(self, st+sname, bool(value))
        # - current status
        self.iterindex = int(f.iterindex)
        self.nitems = int(f.nitems)
        if f.withtime==-1:
           self.withtime  = None
        else:
            self.withtime  = bool(f.withtime)
        if f.withtime:
            self.lasttime = cdtime.s2c(f.lasttime)
        if N.isscalar(f.bin_edges):
            self.bins = None
        else:
            self.bins = N.asarray(f.bin_edges)
            self.nbins = self.bins.shape[0]-1
            self._baxis = f.getAxis('hbin').clone()
        if self.nitems==0: # Still no data
            f.close()
            return 0
        # - already had some data
        self.dual = bool(f.dual)
        self.ns = int(f.ns)
        self.nt = int(f.nt)
        self._nts = f.nts.tolist()
        self.tstats = bool(f.tstats)
        self.sstats = bool(f.sstats)
        if not self.withtime:
            self._stimes = None

        # Spatial statistics
        if self.sstats:

            # Time axes
            if self.withtime:
                self._stimes = tuple([[] for i in xrange(self.nitems)])
                for i, tt in enumerate(self._stimes):
                    taxis = f.getAxis('t'+str(i))
                    tvalues = self._aslist_(taxis[:])
                    oldid = taxis.stataccum_oldid
                    for tvals in tvalues:
                        tx = create_time(tvals, taxis.units, id=oldid)
                        cp_atts(taxis, tx, id=False, exclude=[oldid])
                        self._stimes[i].append(tx)

            # Count
            self._scount = self._load_array_(f, id='scount')

            # Other stats
            self._sstats = {}
            for key in self.single_stats:
                if not self.dual: # single var
                    vid = 's' + key
                    if vid not in f.variables:
                        continue
                    self._sstats[key] = self._load_array_(f, vid),
                else: # two vars
                    for i in xrange(self.nitems):
                        vid = 's%s%s'%(key, str(i))
                        if vid not in f.variables:
                            break
                        self._sstats.setdefault(key, ())
                        self._sstats[key] += self._load_array_(f, vid),
            for key in self.dual_stats:
                vid = 's%s'%key
                if vid in f.variables:
                    self._sstats[key] = self._load_array_(f, vid)

        # Temporal statistics
        if self.tstats:

            # Count
            self._tcount = self._load_array_(f, 'tcount')

            # Other stats
            for key in self._dual_accums+self._single_accums:
                tid = 't'+key
                if not getattr(self, tid): continue
                if key in self._dual_accums:
                    value = self._load_array_(f, tid)
                    setattr(self, '_'+tid, value)
                else:
                    value = ()
                    for i in xrange(self.nitems):
                        value += self._load_array_(f, tid+str(i)),
                    setattr(self, '_'+tid, value)

            # Templates
            # - base arrays
            self._tbase = N.zeros(self.ns)
            if self.thist:
                self._thbase = N.zeros((self.nbins, self.ns), 'l')
            # - cdat templates
            self._ttemplates = ()
            self._thtemplates = None
            if self.thist:
                self._thtemplates = ()
            for i in xrange(self.nitems):
                prefix = 'var%i_'%i
                for vname in f.variables:
                    if vname.startswith(prefix) and vname != prefix+'atts': break
                ttpl = f(vname)
                _rm_id_prefix_(ttpl, 'var%i_'%i, exc=self._baxis)
                self._ttemplates += ttpl,
                if self.thist:
                    self._thtemplates += self._template_t2ht_(ttpl),

        # Attributes
        self._atts = ()
        for ivar in xrange(self.nitems):
            attrs = f['var%i_atts'%ivar].attributes.copy()
            attrs['id'] = attrs['stataccum_id']
            del attrs['stataccum_id']
            self._atts += attrs,

        f.close()
        return self.iterindex




    @staticmethod
    def _base2stats_(size, count, sum, sqr, prod, vmin, vmax, hist,
            domean, dostd, dobias, dorms, docrms, docorr, doavail,
            docount, domin, domax,
            dohist):

        # Denominators
        if isinstance(count, N.ndarray):
            bad = count==0
            count_ = N.where(bad, 1, count)
            countd =  N.where(bad, 0., 1./count_)
#            bad = count<2
#            count = N.where(bad, 2, count)
#            countdm1 =  N.where(bad, 0., 1./(count-1))
#            del count
        else:
            countd = 0. if count == 0 else 1./count
#            countdm1 = 0. if scount < 2 else 1./(count-1)

        # Check dual inputs
        for vname in ('count',' sum',' sqr',' prod',' vmin',' vmax',' hist'):
            val = eval(vname)
            if val is None: continue
            dual = len(val)==2
            if dual: break
        if dual:
            if sum is not None: sum0, sum1 = sum
            if sqr is not None: sqr0, sqr1 = sqr
        else:
            if sum is not None: sum0 = sum[0]
            if sqr is not None: sqr0 = sqr[0]

        # Init
        results = {}

        # Count
        if docount:
            results['count'] = count

        # Availability
        if doavail:
            results['avail'] = 1.* count / size

       # Mean
        if domean:
            results['mean'] = sum0 * countd
            if dual:
                results['mean'] = results['mean'], sum1 * countd


        # Standard deviation & sum(Xi'^2) & sum(Yi'^2)
        if dostd or docorr:
            sqrp0 = - sum0**2 * countd # sum(Xi'^2)
            sqrp0 += sqr0
            if dostd:
                try:
                    results['std'] = N.sqrt(sqrp0 * countd)#m1)
                except:
                    pass
            if not docorr: del sqrp0
            if dual:
                sqrp1 = - sum1**2 * countd # sum(Yi'^2)
                sqrp1 += sqr1
                if dostd:
                    results['std'] = results['std'], N.sqrt(sqrp1 * countd)
                if not docorr: del sqrp1

        # Bias
        if dobias:
            results['bias'] = (sum1-sum0) * countd

        # Absolute RMS
        if dorms:
            rms = sqr0 + sqr1
            rms -= prod * 2
            results['rms'] = N.sqrt(rms * countd)#m1)
            del sqr0, sqr1, rms

        # sum(Xi'Yi')
        if docrms or docorr:
            prodp = prod - sum0 * sum1 * countd

        # Centered RMS
        if docrms:
            crms = sqrp0 + sqrp1
            crms -= prodp * 2
            if not docorr: del prodp
            results['crms'] = N.sqrt(crms * countd)#m1)
            del crms

        # Correlation
        if docorr:
            results['corr'] = prodp / N.sqrt(sqrp0*sqrp1)
            del prodp, sqrp0, sqrp1

        # Min/max
        if domin:
            results['min'] = vmin if dual else vmin[0]
        if domax:
            results['max'] = vmax if dual else vmax[0]

        # Histograms
        if dohist:
            results['hist'] = hist if dual else hist[0]

        return results

    def _checkstat_(self, name):
        if self.iterindex==0:
            raise StatAccumError("Can't get statistics '%s' since still nothing has been accumulated"%name)
        if not getattr(self, name):
            raise StatAccumError("Can't get statistics '%s' since it set to OFF")

    def _get_masks_(self, stats, count):
        masks = {}
        count = N.asarray(count)
        mask0 = count<1
#        mask1 = count<2
        del count
        for key in stats:
            if key in ['mean', 'bias']:
                masks[key] = mask0
            elif key in ['avail', 'count']:
                masks[key] = None
            else:
                masks[key] = mask0
        return masks



    @staticmethod
    def _format_var_(vv, name, mask=None, prefix=True, suffix=True, templates=None,
        htemplates=None, attributes=None, single=True, **kwargs):

        # Some inits
        ist = name.startswith('t')
        name = name[1:]
        if prefix is True:
            prefix = 'Temporal ' if ist else 'Spatial '
        elif not prefix:
            prefix = ''
        kwargs['copy'] = 0

        # Aways have a list
        if isinstance(vv, tuple):
            vv = list(vv)
            single = False
        elif not isinstance(vv, list):
            vv = [vv]
            single = True
        elif single:
            single = len(vv)==1
        else:
            single = False
        dual = len(vv)==2
        if templates is not None and not isinstance(templates, (list, tuple)):
            templates = [templates]*len(vv)
        if htemplates is not None and not isinstance(htemplates, (list, tuple)):
            htemplates = [htemplates]*len(vv)
        if attributes is None:
            attributes = ({}, {})
        elif not isinstance(attributes, (list, tuple)):
            attributes = [attributes]*len(vv)

        # Long_name
        long_name = prefix
        if name=='avail':
            long_name += 'availability'
        if name=='count':
            long_name += 'count'
        elif name=='mean':
            long_name += 'average'
        elif name=='std':
            long_name += 'standard deviation'
        elif name=='bias':
            long_name += 'bias'
        elif name=='rms':
            long_name += 'absolute RMS difference'
        elif name=='crms':
            long_name += 'centered RMS difference'
        elif name=='corr':
            long_name += 'correlation'
        elif name=='hist':
            long_name += 'histogram'

        # Type change
        if 'count' in name:
            dtype = 'l'
        else:
            dtype = None

        for i, v in enumerate(vv):

            # From template or not
            hvar = v.ndim==2 # Histogram-like variables: first axis = bins
            hist = htemplates and hvar
            if templates is not None:

                # From axis (sstats)
                if isaxis(templates[i]):
#                    if hvar and not ist:
#                        v = N.ma.transpose(v) # time first, not bins
                    var = MV2.asarray(v)
                    var.setAxis(0, templates[i])
#                    try:
#                        var.setAxis(int(hist), templates[i])
#                    except:
#                        pass
                    if hist:
                        var.setAxis(1, htemplates[i])

                 # From variable (tstats)
                else:
                    try:
                        var = (htemplates if hvar else templates)[i].clone()
                    except:
                        pass
                    try:
                        var[:] = v.reshape(var.shape)
                    except:
                        pass
            else:
                var = MV2.asarray(v)
            del v

            # Type
            if dtype:
                var = var.astype(dtype)

            # Mask
            if mask is not None:
                func = N.resize if var.size!=mask.size else N.reshape
                var[:] = MV2.masked_where(func(mask, var.shape), var, copy=0)

            # Attributes
            # - id and long_name
            var.id = attributes[i]['id']+'_'+name
            var.long_name = long_name
            if suffix and isinstance(suffix, basestring):
                var.long_name += suffix
            # - single stats
            if name in ['mean', 'std', 'hist', 'min', 'max']:
                if "units" in attributes[i]:
                    var.units = attributes[i]['units']
                if suffix is True and "long_name" in attributes[i]:
                    var.long_name += ' of '+attributes[i]['long_name']
            # - dual stats
            elif name in ['bias', 'rms', 'crms']:
                for attr in attributes: # loop on both variables attributes
                    if "units" in attr and not hasattr(var, 'units'):
                        var.units = attr['units']
                    if suffix is True and "long_name" in attr and var.long_name==long_name:
                        var.long_name += ' of '+attr['long_name']
            vv[i] = var

        if single:
            return vv[0]
        return tuple(vv)


    def get_tstats(self, **kwargs):
        """Get all temporal statistics as a dictionary"""
        if not self.tstats: return {}

        # From cache
        if getattr(self, '_cache_titerindex', -1)==self.iterindex:
            return self._cache_tstats

        # Compute stats
        tstats = self._base2stats_(self.nt, self._tcount, self._tsum if self.tsum else None,
            self._tsqr if self.tsqr else None, self._tprod if self.tprod else None,
            self._tmin if self.tmin else None, self._tmax if self.tmax else None,
            self._thist if self.thist else None,
            self.tmean, self.tstd, self.tbias,
            self.trms, self.tcrms, self.tcorr, self.tavail, self.tcount,
            self.tmin, self.tmax, self.thist)

        # Masks
        masks = self._get_masks_(tstats, self._tcount)

        # Format and cache
        kwargs.update(templates=self._ttemplates, htemplates=self._thtemplates,
            attributes=self._atts, single=True, baxis=self._baxis)
        self._cache_tstats = dict([(name,
            self._format_var_(var, 't'+name, mask=masks[name], **kwargs))
            for name, var in tstats.iteritems()])
        self._cache_titerindex = self.iterindex
        del masks
        return self._cache_tstats

    def get_tavail(self, **kwargs):
        """Get the temporal availability"""
        self._checkstat_('tavail')
        return self.get_tstats()['avail']

    def get_tcount(self, **kwargs):
        """Get the temporal count"""
        self._checkstat_('tcount')
        return self.get_tstats()['count']

    def get_tmean(self, **kwargs):
        """Get the temporal standard deviation"""
        self._checkstat_('tmean')
        return self.get_tstats()['mean']

    def get_tstd(self, **kwargs):
        """Get the temporal mean"""
        self._checkstat_('tstd')
        return self.get_tstats()['std']

    def get_tbias(self, **kwargs):
        """Get the temporal bias"""
        self._checkstat_('tbias')
        return self.get_tstats()['bias']

    def get_trms(self, **kwargs):
        """Get the temporal absolute RMS difference"""
        self._checkstat_('trms')
        return self.get_tstats()['rms']

    def get_tcrms(self, **kwargs):
        """Get the temporal centered RMS difference"""
        self._checkstat_('tcrms')
        return self.get_tstats()['crms']

    def get_tcorr(self, **kwargs):
        """Get the temporal correlation"""
        self._checkstat_('tcorr')
        return self.get_tstats()['corr']

    def get_thist(self, **kwargs):
        """Get the temporal histogram"""
        self._checkstat_('thist')
        return self.get_tstats()['hist']

    def get_sstats(self, **kwargs):
        """Get all spatial statistics as a dictionary"""
        if not self.sstats: return {}

        # From cache
        if getattr(self, '_cache_siterindex', None)==self.iterindex:
            return self._cache_sstats

        # Aggregate pure numeric stats
        sstats = {}
        for stype, vals in self._sstats.items():
            if isinstance(vals, tuple):
                sstats[stype] = self._asarray_(vals[0]), self._asarray_(vals[1])
            else:
                sstats[stype] = self._asarray_(vals)

#        # Availability
#        if self.savail:
#            if self.smask is not None:
#                if isinstance(self.smask , N.ndarray):
#                    maxcount =  (~self.smask).astype('l').sum()
#                else:
#                    maxcount = self.smask
#            else:
#                maxcount = self._scount.max()
#            sstats['avail'] = self._scount*1./maxcount

        # Masks
        masks = self._get_masks_(sstats, self._asarray_(self._scount))

        # Time axes
        stimes = None if self._stimes is None else \
            [MV2_axisConcatenate(stime) for stime in self._stimes]

        # Format and cache
        kwargs.update(templates=stimes, htemplates=self._baxis, attributes=self._atts)
        self._cache_sstats = dict([(name,
            self._format_var_(var, 's'+name, masks[name], **kwargs))
            for name, var in sstats.iteritems()])
        self._cache_siterindex = self.iterindex
        del masks
        return self._cache_sstats

    def get_savail(self, **kwargs):
        """Get the spatial availability"""
        self._checkstat_('savail')
        return self.get_sstats()['avail']

    def get_scount(self, **kwargs):
        """Get the spatial count"""
        self._checkstat_('scount')
        return self.get_sstats()['count']

    def get_smean(self, **kwargs):
        """Get the spatial standard deviation"""
        self._checkstat_('smean')
        return self.get_sstats()['mean']

    def get_sstd(self, **kwargs):
        """Get the spatial mean"""
        self._checkstat_('sstd')
        return self.get_sstats()['std']

    def get_sbias(self, **kwargs):
        """Get the spatial bias"""
        self._checkstat_('sbias')
        return self.get_sstats()['bias']

    def get_srms(self, **kwargs):
        """Get the spatial absolute RMS difference"""
        self._checkstat_('srms')
        return self.get_sstats()['rms']

    def get_scrms(self, **kwargs):
        """Get the spatial centered RMS difference"""
        self._checkstat_('scrms')
        return self.get_sstats()['crms']

    def get_scorr(self, **kwargs):
        """Get the spatial correlation"""
        self._checkstat_('scorr')
        return self.get_sstats()['corr']

    def get_shist(self, **kwargs):
        """Get the spatial histogram"""
        self._checkstat_('shist')
        return self.get_sstats()['hist']

    def get_stats(self):
        """Get all statistics as dict"""
        stats = {}
        for key, val in self.get_tstats().items():
            stats['t'+key] = val
        for key, val in self.get_sstats().items():
            stats['s'+key] = val
        return stats

    def __getitem__(self, key):
        return self.get_stats()[key]

def _get_var_and_axes_(var, exc=None):
    """Get a list of a variable and its axes"""
    if exc is None: exc = []
    if not isinstance(exc, (list, tuple)):
        exc = [exc]
    out = [var]+var.getAxisList()
    if var.getGrid() and var.getLongitude().id!=var.getAxis(-1).id:
        out.extend([var.getLongitude(), var.getLatitude()])
    if not exc: return out
    exc0 = map(id, exc)
    exc1 = [e.id for e in exc]
    return filter(lambda o: id(o) not in exc0 and o.id not in exc1, out)

def _rm_id_prefix_(var, prefix, exc=None):
    """Remove prefix in ids of a variable and its axes"""
    for obj in _get_var_and_axes_(var, exc=exc):
        if obj.id.startswith(prefix):
            obj.id = obj.id[len(prefix):]

def _add_id_prefix_(var, prefix, exc=None):
    """Add prefix in ids of a variable and its axes when not present"""
    for obj in _get_var_and_axes_(var, exc=exc):
        if not obj.id.startswith(prefix):
            obj.id = prefix+obj.id
