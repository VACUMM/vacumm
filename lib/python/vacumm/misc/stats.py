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

from genutil.salstat import betai,_gammaln
from genutil.statistics import percentiles
import numpy as N,numpy.ma as MA,MV2
from vacumm import VACUMMError
from vacumm.misc.axes import isaxis
from vacumm.misc import MV2_axisConcatenate

__all__ = ['corr_proba', 'ensrank', 'qtmin', 'qtmax', 'qtminmax', 'StatAccum', 'StatAccumError']
__all__.sort()

def corr_proba(r, ndata, ndataset=2, dof=False):
    """Probability of rejecting correlations

    - **r**: Correlation coefficient
    - **ndata**: Number of records use for correlations
    - **ndataset**, optional:  Number of datasets (1 for autocorrelations, else 2) [default: 2]

    .. todo::

        This must be rewritten using :mod:`scipy.stats`
    """

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

    _single_stats = 'avail', 'mean', 'std'
    _dual_stats = 'bias', 'rms', 'crms', 'corr'
    _all_stats = _single_stats + _dual_stats
    def __init__(self,
        tall=False, tavail=None, tmean=None, tstd=None, tbias=None, trms=None, tcrms=None, tcorr=None,
        sall=False, savail=None, smean=None, sstd=None, sbias=None, srms=None, scrms=None, scorr=None,
#        tmask=None, smask=None,
        withtime=None):
        for st in 'st':
            activate = False
            for stype in self._all_stats:
                value = eval(st+stype)
                if value is None: value = eval(st+'all')
                setattr(self, st+stype, value)
                activate |= value
            setattr(self, st+'stats', activate)
#        self.tmask = tmask
#        self.smask = smask
        self.iterindex = 0
        self.withtime = withtime

    def append(self, *items):
        """Append data to the accumulator"""
        # Initialization
        if self.iterindex==0:

            # Two variables?
            self.dual = len(items)>1

            # Do we have time in input variables
            if self.withtime is None:
                self.withtime = items[0].getTime() is not None

            # Init template for output temporal statistics
            if self.tstats:
                self._ttemplates = []
                selspace = 0 if self.withtime else slice(None)
                for item in items:
                    tpl = item[selspace].clone()
                    for att in tpl.attributes:
                        delattr(tpl, att)
                    self._ttemplates.append(tpl)
                self._tbase = items[0].filled()[selspace]*0.
                self._tbase.shape = self._tbase.size,

            # Check spatial statistics
            if self.sstats and items[0].ndim<2:
                self.sstats = self.smean = self.sstd = self.sbias = self.srms = self.scrm = self.scorr = self.savail = False

            # Save some attributes
            self._atts = []
            for var in items:
                attr = {}
                for att in 'id', 'units', 'long_name', 'standard_name':
                    if hasattr(var, att):
                        attr[att] = getattr(var, att)
                self._atts.append(attr)


            # Flags
            if not self.dual:
                self.tbias = self.trms = self.tcrms = self.tcorr = False
                self.sbias = self.srms = self.scrms = self.scorr = False
            self.tstats = self.tmean or self.tstd or self.tbias or self.trms or self.tcrms or self.tavail
            self.sstats = self.smean or self.sstd or self.sbias or self.srms or self.scrms or self.savail
            self.tsum = self.tmean or self.tbias or self.tstd or self.tcrms or self.tcorr
            self.tsqr = self.tmean or self.tstd or self.trms or self.tcrms or self.tcorr
            self.tprod = self.trms or self.tcrms or self.tcorr
            self.ssum = self.smean or self.sbias or self.sstd or self.scrms or self.scorr
            self.ssqr = self.smean or self.sstd or self.srms or self.scrms or self.scorr
            self.sprod = self.srms or self.scrms or self.scorr

            # Init accumulated variables

            # - temporal statistics
            if self.tstats:
                if self.tsum:
                    self._tsum = [self._tbase.copy()]
                    if self.dual: self._tsum.append(self._tbase.copy())
                # - square
                if self.tsqr:
                    self._tsqr = [self._tbase.copy()]
                    if self.dual: self._tsqr.append(self._tbase.copy())
                # - product
                if self.tprod: self._tprod = self._tbase.copy()
                # - count
                self._tcount = self._tbase.copy().astype('l')
                self._tsize = 0l

            # - spatial statistics
            if self.sstats:
                self._sstats = {}
#                for name in self._all_stats:
#                    if not getattr(self, 's'+name): continue
#                    self._sstats[name] = []
#                    setattr(self, '_s'+name, [])
                self._stimes = [[] for i in xrange(len(items))] if self.withtime else None
                self._scount = []
                self._ssize = items[0].size/items[0].shape[0]

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
        item0 = items[0].asma().reshape(nt, -1)
        mask0 = N.ma.getmaskarray(item0)
        if self.dual:
            item1 = items[1].asma().reshape(nt, -1)
            mask1 = N.ma.getmaskarray(item1)
            mask = mask0 | mask1
        else:
            mask = mask0
        item0[mask] = N.ma.masked
        item0 = item0.filled(0.)
        if self.dual:
            item1[mask] = N.ma.masked
            item1 = item1.filled(0.)

        # Accumulate

        # - temporal statistics
        if self.tstats:
            self._tcount += (~mask).astype('l').sum(axis=0)
            self._tsize += items[0].shape[0]
            if self.tsum:
                self._tsum[0] += item0.sum(axis=0)
                if self.dual: self._tsum[1] += item1.sum(axis=0)
            if self.tsqr:
                self._tsqr[0] += (item0**2).sum(axis=0)
                if self.dual: self._tsqr[1] += (item1**2).sum(axis=0)
            if self.tprod:
                self._tprod += (item0*item1).sum(axis=0)

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

            # Stats
            sstats = self._base2stats_(self._ssize, scount, ssum if self.ssum else None,
                ssqr if self.ssqr else None, sprod if self.sprod else None,
                self.smean, self.sstd, self.sbias, self.srms, self.scrms, self.scorr, self.savail)
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

        self.iterindex += 1


    def __iadd__(self, items):
        if isinstance(items, tuple):
            self.append(*items)
        else:
            self.append(items)
        return self


    def dump(self, ncfile):
        f = cdms2.open('ncfile', 'w')

        # Scalars
        for sname in _all_stats + ['sum', 'sqr', 'prod', 'stats']:
            for st in 'st':
                setattr(f, st+sname, getattr(self, st+sname))
        f.iterindex = self.iterindex
        f.withtime = self.withtime

        # Arrays
        for st in 'st':
            if getattr(self, st+'stats'):
                count = MV2.array(getattr(self, '_%scount'%st), id=st+'count', copy=0)
                count.setAxis(0, axis)
                f.write(count)
                axis = cdms2.createAxis(N.arange(len(count)), id='s' if st=='t' else 't')
                stats = getattr(self, '_%sstats'%st)
                for key in stats:
                    var = MV2.array(stats[key], copy=0, id=st+key)
                    var.setAxis(0, axis)
                    f.write(var)

    @staticmethod
    def _base2stats_(size, count, sum, sqr, prod,
        domean, dostd, dobias, dorms, docrms, docorr, doavail):

        # Denominators
        if isinstance(count, N.ndarray):
            bad = count==0
            count = N.where(bad, 1, count)
            countd =  N.where(bad, 0., 1./count)
#            bad = count<2
#            count = N.where(bad, 2, count)
#            countdm1 =  N.where(bad, 0., 1./(count-1))
#            del count
        else:
            countd = 0. if count == 0 else 1./count
#            countdm1 = 0. if scount < 2 else 1./(count-1)

        # Dual inputs
        dual = len(sum)==2 or len(sqr)==2
        if dual:
            sum0, sum1 = sum
            sqr0, sqr1 = sqr
        else:
            sum0 = sum[0]
            sqr0 = sqr[0]

        # Init
        results = {}

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
                results['std'] = N.sqrt(sqrp0 * countd)#m1)
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
            elif key in ['avail']:
                masks[key] = None
            else:
                masks[key] = mask0
        return masks



    @staticmethod
    def _format_var_(vv, name, mask=None, prefix=True, suffix=True, templates=None,
        attributes=None, single=True, **kwargs):
        if name=='tstd':
            pass
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
        if attributes is not None and not isinstance(attributes, (list, tuple)):
            attributes = [attributes]*len(vv)

        # Long_names
        long_name = prefix
        if name=='avail':
            long_name += 'availability'
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

        for i, v in enumerate(vv):

            # From template or not
            if templates is not None:

                # From axis (sstats)
                if isaxis(templates[i]):
                    var = MV2.asarray(v)
                    try:
                        var.setAxis(0, templates[i])
                    except:
                        pass
                 # From variable (tstats)
                else:
                    var = templates[i].clone()
                    var[:] = v.reshape(var.shape)
            else:
                var = MV2.asarray(v)
            del v

            # Mask
            if mask is not None:
                mask.shape = var.shape
                var[:] = MV2.masked_where(mask, var, copy=0)

            # Attributes
            if attributes is None:
                attr = {}
            else:
                attr = attributes[i]
            var.id = attr['id']+'_'+name
            var.long_name = long_name
            if suffix and isinstance(suffix, basestring):
                var.long_name += suffix
            if name in ['mean', 'std']:
                if "units" in attr:
                    var.units = attr['units']
                if suffix is True and "long_name" in attr:
                    var.long_name += ' of '+attr['long_name']
            elif name in ['bias', 'rms', 'crms']:
                for att in attr:
                    if "units" in att and not hasattr(var, 'units'):
                        var.units = att['units']
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
        if getattr(self, '_cache_titerindex', None)==self.iterindex:
            return self._cache_tstats

        # Compute stats
        tstats = self._base2stats_(self._tsize, self._tcount, self._tsum if self.tsum else None,
            self._tsqr if self.tsqr else None, self._tprod if self.tprod else None,
            self.tmean, self.tstd, self.tbias,
            self.trms, self.tcrms, self.tcorr, self.tavail)

#        # Availability
#        if self.tavail:
##            if self.tmask is not None:
##                if isinstance(self.tmask , N.ndarray):
##                    maxcount =  (~self.tmask).astype('l').sum()
##                else:
##                    maxcount = self.tmask
##            else:
##                maxcount = self._tcount.max()
#            tstats['avail'] = self._tcount*1./self._tsize

        # Masks
        masks = self._get_masks_(tstats, self._tcount)

        # Format and cache
        kwargs.update(templates=self._ttemplates, attributes=self._atts, single=True)
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

    def get_sstats(self, **kwargs):
        """Get all spatial statistics as a dictionary"""
        if not self.sstats: return {}

        # From cache
        if getattr(self, '_cache_siterindex', None)==self.iterindex:
            return self._cache_sstats

        # Compute stats
        sstats = {}
        for stype, vals in self._sstats.items():
            if isinstance(vals, tuple):
                sstats[stype] = N.concatenate(vals[0]), N.concatenate(vals[1])
            else:
                sstats[stype] = N.concatenate(vals)

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
        masks = self._get_masks_(sstats, self._scount)

        # Time axes
        stimes = None if self._stimes is None else \
            [MV2_axisConcatenate(stime) for stime in self._stimes]

        # Format and cache
        kwargs.update(templates=stimes, attributes=self._atts)
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

    def get_stats(self):
        """Get all statistics as dict"""
        stats = {}
        for key, val in self.get_tstats().items():
            stats['t'+key] = val
        for key, val in self.get_sstats().items():
            stats['s'+key] = val
        return stats
