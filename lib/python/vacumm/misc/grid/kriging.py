# -*- coding: utf8 -*-
"""
Kriging utilities inspired from the AMBHAS library (http://www.ambhas.com/).
"""
# Copyright or © or Copr. Actimar/IFREMER (2013-2016)
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

import gc
import multiprocessing
from multiprocessing import Pool,  cpu_count
import warnings

import numpy as N
import pylab as P

if not hasattr(N, 'isclose'):
    from vacumm.misc import closeto as isclose
else:
    isclose = N.isclose
from vacumm.misc import bivariate_normal

#from scipy.optimize import curve_fit
def get_blas_func(name):
    try:
        import scipy.linalg.blas
        func = scipy.linalg.blas.get_blas_funcs(name)
    except:
        import scipy.linalg.fblas
        func = getattr(scipy.linalg.fblas, 'd'+name)
    return func
blas_dgemv = get_blas_func('gemv')
def dgemv(a, x): return blas_dgemv(1., a, x)
try:
    from _blaslapack import symm, sytri
except:
#    print 'Falling back to builtin functions'
    dgemm = get_blas_func('gemm')
    def symm(a, b): return dgemm(1., a, b)
    sytri = N.linalg.pinv

from ...misc.misc import kwfilter
from .misc import get_distances

class KrigingError(Exception):
    pass

#: Variogram model types
VARIOGRAM_MODEL_TYPES = ['linear', 'exponential', 'spherical', 'gaussian']
VARIOGRAM_MODEL_TYPES = VARIOGRAM_MODEL_TYPES

#: Default variogram model type
DEFAULT_VARIOGRAM_MODEL_TYPE = 'exponential'

def variogram_model_type(mtype=None):
    """Check the the variogram model type

    :Params:

        - **mtype**, optional: ``None``, and index or a string matching
          an element of :data:`VARIOGRAM_MODEL_TYPES`.
          If set to ``None``, it defaults to :data:`DEFAULT_VARIOGRAM_MODEL_TYPE`.
    """
    if mtype is True:
        return VARIOGRAM_MODEL_TYPES
    errmsg = []
    for i, vtype in enumerate(VARIOGRAM_MODEL_TYPES):
        errmsg += '"%s" (=%i)'%(vtype, i)
    errmsg = 'Invalid variogram model type. Please choose one of: '+ ', '.join(errmsg)
    if mtype is None:
        mtype = DEFAULT_VARIOGRAM_MODEL_TYPE
    if isinstance(mtype, int):
        if i<0 or i>len(VARIOGRAM_MODEL_TYPES)-1:
            raise KrigingError(errmsg)
        return VARIOGRAM_MODEL_TYPES[i]
    if not isinstance(mtype, basestring):
        raise KrigingError(errmsg)
    for vtype in VARIOGRAM_MODEL_TYPES:
        if vtype.startswith(mtype): return vtype
    raise KrigingError(errmsg)


def variogram_model(mtype, n, s, r, nrelmax=0.2):
    """Get the variogram model function from its name"""

    mtype = variogram_model_type(mtype)

    n = max(n, 0)
    n = min(n, nrelmax*s)
    r = max(r, 0)
    s = max(s, 0)

    if mtype == 'linear':
        return lambda h: n + (s-n) * ((h/r)*(h<=r) + 1*(h>r))

    if mtype=='exponential':
        return lambda h: n + (s-n) * (1 - N.exp(-3*h/r))

    if mtype=='spherical':
        return lambda h: n + (s-n)*((1.5*h/r - 0.5*(h/r)**3)*(h<=r) + 1*(h>r))

    if mtype=='gaussian':
        return lambda h: n + (s-n)*(1-N.exp(-3*h**2/r**2))

    raise KrigingError(errmsg)

class VariogramModel(object):
    """Class used when fitting a variogram model to data to better control params"""
    param_names = list(variogram_model.func_code.co_varnames[1:])
    param_names.remove('nrelmax')
    def __init__(self, mtype, **kwargs):
        self.mtype = variogram_model_type(mtype)
        self.fixed_params = dict([(p, v) for (p, v) in kwargs.items()
            if p in self.param_names and v is not None])

    def get_all_kwargs(self, pp):
        """Get arguments list to :func:`variogram_model` by merging variable params `p`
        and :attr:`fixed_params`
        """
        pp = list(pp)
        return dict([(p, (self.fixed_params[p] if p in self.fixed_params else pp.pop(0)))
            for p in self.param_names])

    def get_var_args(self, **kwargs):
        """Get variable arguments list from specified params

        .. note:: Result cannot contain ``None``
        """
        vargs = [kwargs[p] for p in self.param_names if p not in self.fixed_params]
        if None in vargs:
            raise VariogramModelError('Variable arguments cannot contains Nones')
        return vargs

    def __call__(self, d, *pp):
        """Call the variogram model function"""
        return self.get_variogram_model(pp)(d)

    def get_variogram_model(self, pp):
        """Get the variogram model function using `pp` variable arguments"""
        kwargs = self.get_all_kwargs(pp)
        return variogram_model(self.mtype, **kwargs)


class VariogramModelError(KrigingError):
    pass

def _get_xyz_(x, y, z=None, check=True, noextra=True, getmask=False):
    if not x.ndim==1 and x.shape!=y.shape :
        raise KrigingError('x, y must have the same 1D shape')
    if z is not None:
        if (noextra or z.ndim==1) and x.shape!=z.shape:
                raise KrigingError('z must have the same 1D shape as x and y')
        elif not noextra and z.ndim!=1 and (z.ndim!=2 or x.shape!=z.shape[1:]):
            raise KrigingError('z must have 2 dims with its last dim of same length as x and y')
    mask = N.ma.nomask
    if N.ma.isMA(z): mask = N.ma.getmaskarray(z)
    if N.ma.isMA(x): mask |= N.ma.getmaskarray(x)
    if N.ma.isMA(y): mask |= N.ma.getmaskarray(y)
    if mask is not N.ma.nomask and mask.any():
        if mask.shape==2: # permanent mask
            mask = N.logical_and.reduce(mask, axis=0)
        good = ~mask
        if check and not good.any():
            raise KrigingError('All data are masked')
        try:
            x = x.compress(good)
        except:
            pass
        y = y.compress(good)
        z = z.compress(N.resize(good, z.shape))
    res = x, y, z
    if getmask:
        res += mask,
    return res

def variogram(x, y, z, binned=None, nmax=1500, nbindef=30, nbin0=None,
        nbmin=10, distmax=None,  distfunc='simple', errfunc=None):
    """Estimate variogram from data

    :Params:

        - **x/y/z**: 1D arrays of positions and data.
        - **nmax**, optional: Above this number, size of
          the sampe is reduced using undersampling.
        - **binned**, optional: If set to a number,
          data are arranged in bins to estimate
          variogram. If set to ``None``, data are
          arranged in bins if the number of pairs
          of points is greater than ``nbindef*nbmin``.
        - **nbindef**, optional: Default number
          of bins (not used if ``binned`` is a number).
        - **nbin0**, optional: If set to a number > 1,
          the first bin is split into nbin0 sub-bins.
          If set to ``None``, it is evaluated with
          ``min(bins[1]/nbmin, nbin)``.
        - **nbmin**, optional: Minimal number of points
          in a bin.
        - **distmax**, optional: Max distance to consider.
        - **distfunc**: Function to compute distances, or a mode argument to
          :func:`~vacumm.misc.grid.misc.get_distances`.
        - **errfunc**, optional: Callable function to compute "errors" like square
          root difference between to z values. It take two arguments and
          defaults to :math:`(z1-z0)^2/2`.

    """
    x, y, z = _get_xyz_(x, y, z)
    npts = x.shape[0]

    # Undepsample?
    if npts>2*nmax:
        samp = npts/nmax
        x = x[::samp]
        y = y[::samp]
        z = z[::samp]
        npts = x.shape[0]

    # Distances
    dd = get_distances(x, y, x, y, mode=distfunc)

    # Variogram
    if errfunc is None:
        errfunc = lambda a0, a1: 0.5*(a1-a0)**2
    vv = errfunc(*N.meshgrid(z, z))

    # Unique
    ii = N.indices(dd.shape)
    iup = ii[1]>ii[0]
    d = dd[iup]
    v = vv[iup]
    del dd, vv

    # Max distance
    if distmax:
        valid = d<=distmax
        d = d[valid]
        v = v[valid]
        del valid

    # Direct variogram?
    if binned is None and len(d)>nbindef*nbmin: binned = True
    if binned is True: binned = nbindef
    if not binned: return d, v

    # Rebin
    # - first try
    ii = N.argsort(d)
    np = d.shape[0]
    nbin = binned
    edges = N.linspace(0, np-1, nbin+1).astype('i').tolist()
    # - more details in first bin
    if nbin0 is None: # do we need more details?
        NBIN0MAX = 10
        nbin0 = min(edges[1]/nbmin, NBIN0MAX)
    if nbin0>1: # split first bin
        edges = N.linspace(0., edges[1], nbin0+1).astype('i')[:-1].tolist()+edges[1:]
        nbin = nbin - 1 + nbin0 # len(edges)-1
    # - rebinning
    db = N.empty(nbin)
    vb = N.empty(nbin)
    for ib in xrange(nbin):
        iib = ii[edges[ib]:edges[ib+1]+1]
        db[ib] = d[iib].mean()
        vb[ib] = v[iib].mean()
    return db, vb

def variogram_fit(x, y, z, mtype=None, getall=False, getp=False, geterr=False,
        distfunc='simple', errfunc=None, **kwargs):
    """Fit a variogram model to data and return the function

    :Example:

        >>> vm, errs = variogram_fit(x, y, z, 'linear', n=0, distmax=30e3, geterr=True)

    :Params:

        - **x/y/z**: Position and data.
        - **mtype**: Variogram model type (see ::`variogram_model_type`).
        - **getall**: Get verything in a dictionary whose keys are

            - ``"func"``: model function,
            - ``"err"``: fitting error,
            - ``"params"``: all parameters has a dictionary,
            - ``"popt"``: parameters than where optimised,
            - ``vm"``: :class:`VariogramModel` instance,
            - ``"mtype"``: variogram model type.

        - **getp**, optional: Only return model parameters. Return them as
          a `class:`dict` if equal to ``2``.
        - **variogram_<param>**, optional: ``param`` is passed to :func:`variogram`.
        - **distfunc**: Function to compute distances, or a mode argument to
          :func:`~vacumm.misc.grid.misc.get_distances`.
        - **errfunc**, optional: Callable function to compute "errors" like square
          root difference between to z values. It take two arguments and
          defaults to :math:`\sqrt(z1^2-z0^2)/2`.

          .. warning:: use "haversine" if input coordinates are in degrees.

        - Extra keywords are those of :func:`variogram_model`.
          They can be used to freeze some of the parameters.

          >>> variogram_fit(x, y, z, mtype, n=0) # fix the nugget

    """
    kwv = kwfilter(kwargs, 'variogram_')
    kwv.setdefault("distfunc", distfunc)
    kwv.setdefault("errfunc", errfunc)

    # Estimated variogram
    d, v = variogram(x, y, z, **kwv)

    # Variogram model
    vm = VariogramModel(mtype, **kwargs)

    # First guess of paramaters
    imax = N.ma.argmax(v)
    p0 = vm.get_var_args(n=0., s=v[imax], r=d[imax])

    # Fitting
#    p, e = curve_fit(vm, d, v, p0=p0) # old way: no constraint
    from scipy.optimize import minimize
    func = lambda pp: ((v-vm.get_variogram_model(pp)(d))**2).sum()
    warnings.filterwarnings('ignore', 'divide by zero encountered in divide')
    p = minimize(func, p0, bounds=[(N.finfo('d').eps, None)]*len(p0),
        method='L-BFGS-B')['x']
    del warnings.filters[0]

    # Output
    if getall:
        return dict(
            func=vm.get_variogram_model(p),
            err=(vm.get_variogram_model(p)(d)-v).std(),
            params = vm.get_all_kwargs(p),
            popt = p)
    if int(getp)==2:
        res = vm.get_all_kwargs(p)
    elif getp:
        res = p
    else:
        res = vm.get_variogram_model(p)
    if not geterr:
        return res
    return res,  (vm.get_variogram_model(p)(d)-v).std()

def variogram_multifit(xx, yy, zz, mtype=None, getall=False, getp=False, **kwargs):
    """Same as :func:`variogram_fit` but with several samples"""
    vm = VariogramModel(mtype, **kwargs)
    pp = []
    for i, (x, y, z) in enumerate(zip(xx, yy, zz)):
        x, y, z = _get_xyz_(x, y, z, check=False)
        if len(x)==0: continue
        p = variogram_fit(x, y, z, mtype, getp=True, **kwargs)
        pp.append(p)
    pp = N.asarray(pp)
    if pp.shape[0]==0:
        raise KrigingError('All data are masked')
    mp = N.median(pp, axis=0)
    if getall:
        return dict(
            func=vm.get_variogram_model(mp),
            err=None,
            params = vm.get_all_kwargs(mp),
            popt = mp)
    if int(getp)==2:
        return vm.get_all_kwargs(mp)
    if getp:
        return mp
    return vm.get_variogram_model(mp)

def cloud_split(x, y, npmax=1000, getdist=True, getcent=True):
    """Split data intot cloud of points of max size npmax:

    :Returns: ``None`` if ``len(x)<=npmax``

        Else ``indices``
        or ``(indices, global_distorsion, distortions)``.

    """
    from scipy.cluster.vq import kmeans, vq

    # Nothing to do
    csize = len(x)
    if npmax<2 or csize<=npmax: return

    # Loop on the number of clusters
    nclust = 2
    points = N.vstack((x, y)).T
    ii = N.arange(csize)
    while csize > npmax:
        centroids, global_distorsion = kmeans(points, nclust)
        indices, distorsions = vq(points, centroids)
        sindices = [ii[indices==nc] for nc in xrange(nclust)]
        csizes = [sii.shape[0] for sii in sindices]
        order = N.argsort(csizes)[::-1]
        csize = csizes[order[0]]
        if getdist: sdistorsions = [distorsions[sii] for sii in sindices]
        nclust += 1

    #  Reorder
    sindices = [sindices[i] for i in order]
    if getdist:
        sdistorsions = [sdistorsions[i] for i in order]
        dists = global_distorsion,  sdistorsions
    if getcent:
        centroids = centroids[order]

    # Output
    if not getdist and not getcent: return indices
    ret = sindices,
    if getcent: ret += centroids,
    if getdist: ret += dists,
    return ret

def syminv(A):
    """Invert a symetric matrix

    :Params:

        - **A**: (np+1,np+1) for variogram matrix

    :Return: ``Ainv(np+1,np+1)``

    :Raise: :exc:`KrigingError`

    """
    res = sytri(A.astype('d'))
    if isinstance(res, tuple):
        info = res[1]
        if info: raise KrigingError('Error during call to Lapack DSYTRI (info=%i)'%info)
        return res[0]
    else:
        return res

class CloudKriger(object):
    """Ordinary kriger using mutliclouds of points

    Big input cloud of points (size > ``npmax``)
    are split into smaller
    clouds using cluster analysis of distance with
    function :func:`cloud_split`.

    The problem is solved in this way:

        #. Input points are split in clouds if necessary.
        #. The input variogram matrix is inverted
           for each cloud, possibly using
           :mod:`multiprocessing` if ``nproc>1``.
        #. Value are computed at output positions
           using each the inverted matrix of cloud.
        #. Final value is a weighted average of
           the values estimated using each cloud.
           Weights are inversely proportional to the inverse
           of the squared error.

    :Params:

        - **x/y/z**: Input positions and data (masked array).
        - **mtype**, optional: Variogram model type (defaults to 'exp').
          See :func:`variogram_model_type` and :func:`variogram_model_type`.
        - **vgf**, optional: Variogram function. If not set,
          it is estimated using :meth:`variogram_fit`.
        - **npmax**, optional: Maxima size of cloud.
        - **nproc**, optional: Number of processes to use
          to invert matrices. Set it to a number <2 to switch off
          parallelisation.
        - **exact**, optional: If True, variogram is exactly zero when distance is zero.
        - **distfunc**: Function to compute distances, or a mode argument to
          :func:`~vacumm.misc.grid.misc.get_distances`.
        - **errfunc**, optional: Callable function to compute "errors" like square
          root difference between to z values. It take two arguments and
          defaults to :math:`\sqrt(z1^2-z0^2)/2`.
        - Extra keywords are  parameters to the :func:`variogram_model` that must not be
          optimized by :func:`variogram_model`. For instance ``n=0`` fix the
        - Extra keywords are the parameters to the :func:`variogram_model` that must not be
          optimized by :func:`variogram_model`. For instance ``n=0`` fixes the
          nugget to zero.
          This is used only if ``vfg`` is not passed as an argument.


    :Attributes: :attr:`x`, :attr:`y`, :attr:`z`, :attr:`np`,
        :attr:`xc`, :attr:`yc`, :attr:`zc`,  :attr:`npc`,
        :attr:`variogram_function`, :attr:`Ainv`, :attr:`npmax`, :attr:`nproc`.

        .. attribute:: x

            List of all input x positions.

        .. attribute:: y

            List of all input y positions.

        .. attribute:: z

            List of all input data.

        .. attribute:: xc

            List of input x positions of each cloud.

        .. attribute:: yc

            List input of y positions of each cloud.

        .. attribute:: zc

            List of input data of each cloud.
    """

    def __init__(self, x, y, z, krigtype, mtype=None, vgf=None, npmax=1000,
            nproc=None, exact=False, distfunc='simple', errfunc=None,
            mean=None, farvalue=None, **kwargs):
        if krigtype is None:
            krigtype = 'ordinary'
        krigtype = str(krigtype).lower()
        assert krigtype in ['simple', 'ordinary'], ('krigtype must be either '
            '"simple" or "ordinary"')
        self.x, self.y, self.z, self.mask = _get_xyz_(x, y, z, noextra=False, getmask=True)
        self.np = self.x.shape[0]
        self.nt = 0 if self.z.ndim==1 else z.shape[0]
        self.mtype = variogram_model_type(mtype)
        self.npmax = npmax
        if nproc is None:
            nproc = cpu_count()
        else:
            nproc = max(1, min(cpu_count(), nproc))
        self._setup_clouds_()
        self.nproc = min(nproc, self.ncloud)
        if callable(vgf):
            self.variogram_func = vgf
        self._kwargs = kwargs
        self.variogram_fitting_results = None
        self.exact = exact
        self.distfunc = distfunc
        self.errfunc = errfunc
        self.krigtype = krigtype
        self._simple = self.krigtype == 'simple'
        if not self._simple:
            mean = 0.
        elif mean is None:
            mean = self.z.mean()
        self.mean = mean
        self.farvalue = farvalue

    def __len__(self):
        return self.x.shape[0]

    def _get_xyz_(self, x=None, y=None, z=None):
        if x is None: x = self.x
        if y is None: y = self.y
        if z is None: z = self.z
        return _get_xyz_(x, y, z, noextra=False)

    def _setup_clouds_(self):
        """Setup cloud spliting

        Estimate:

            #. The number of procs to use (sef.nproc).
            #. Cloud positions and data (self.xc/yc/zc[ic]).
            #. Weight function for each cloud (self.wc[ic](x,y)).
        """

        # Split?
        self.npc = []
        if self.npmax>2 and self.x.shape[0]>self.npmax:

            # Split in clouds
            indices, centroids, (gdist, dists) = cloud_split(self.x, self.y,
                npmax=self.npmax, getdist=True, getcent=True)
            self.ncloud = len(indices)

            # Loop on clouds
            self.xc = []
            self.yc = []
            self.zc = []
            for ic in xrange(self.ncloud):

                # Positions and data
                for xyz in 'x', 'y', 'z':
                    getattr(self, xyz+'c').append(getattr(self, xyz)[..., indices[ic]].T)

                # Size
                self.npc.append(len(indices[ic]))

        else: # Single cloud

            self.xc = [self.x]
            self.yc = [self.y]
            self.zc = [self.z.T]
            self.ncloud = 1
            self.npc = [self.np]
            self.cwfunc = [lambda x, y: 1.]

    def plot_clouds(self, marker='o', **kwargs):
        """Quickly Plot inputs points splitted in clouds"""
        P.figure()
        for x, y in zip(self.xc, self.yc):
            P.plot(x, y, marker, **kwargs)
        P.show()

    def variogram_fit(self, x=None, y=None, z=None, **kwargs):
        """Estimate the variogram function by using :func:`variogram_fit`"""
        kw = self._kwargs.copy()
        kw.update(kwargs)
        kw['distfunc'] = self.distfunc
        kw['errfunc'] = self.errfunc
        x, y, z = self._get_xyz_(x, y, z)
        if z.ndim==2:
            ne = z.shape[0]
            res = variogram_multifit([x]*ne, [y]*ne, z, self.mtype, getall=True, **kw)
        else:
            res = variogram_fit(x, y, z, self.mtype, getall=True, **kw)
        self.variogram_func = res['func']
        self.variogram_fitting_results = res
        return self.variogram_func

    def get_sill(self):
        vgf = self.variogram_func
        if self.variogram_fitting_results is not None:
            return self.variogram_fitting_results['params']['s']
        return vgf(1e60)

    sill = property(fget=get_sill, doc="Sill")

    def set_variogram_func(self, vgf):
        """Set the variogram function"""
        if not callable(vgf):
            raise KrigingError("Your variogram function must be callable")
        reset = getattr(self, '_vgf', None) is not vgf
        self._vgf = vgf
        if reset: del self.Ainv
    def get_variogram_func(self):
        """Get the variogram function"""
        if not hasattr(self, '_vgf'):
            self.variogram_fit()
        return self._vgf
    def del_variogram_func(self):
        """Delete the variogram function"""
        if hasattr(self, '_vgf'):
            del self._vgf
    variogram_func = property(get_variogram_func, set_variogram_func,
        del_variogram_func, "Variogram function")


    def get_Ainv(self):
        """Get the inverse of A"""

        # Already computed
        if hasattr(self, '_Ainv'):
            return self._Ainv

        # Variogram function
        vgf = self.variogram_func

        # Loop on clouds
        if not hasattr(self, '_dd'): self._dd = []
        Ainv = []
        AA = []
        next = int(not self._simple)
        for ic in xrange(self.ncloud):

            # Get distance between input points
            if len(self._dd)<ic+1:
                dd = get_distances(self.xc[ic], self.yc[ic],
                    self.xc[ic], self.yc[ic], mode=self.distfunc)
                self._dd.append(dd)
            else:
                dd = self._dd[ic]

            # Form A
            np = self.npc[ic]
            A = N.empty((np+next, np+next))
            A[:np, :np] = vgf(dd)
            if self.exact:
                N.fill_diagonal(A, 0)
                A[:np, :np][isclose(A[:np, :np], 0.)] = 0.
            if not self._simple:
                A[-1] = 1
                A[:, -1] = 1
                A[-1, -1] = 0

            # Invert for single cloud
            if self.nproc==1:
                Ainv.append(syminv(A))
            else:
                AA.append(A)

        # Multiprocessing inversion
        if self.nproc>1:
            pool = Pool(self.nproc)
            Ainv = pool.map(syminv, AA, chunksize=1)
            pool.close()

        # Fortran arrays
        Ainv = [N.asfortranarray(ainv, 'd') for ainv in Ainv]
        self.Ainv = Ainv
        return Ainv

    def set_Ainv(self, Ainv):
        """Set the invert of A"""
        self._Ainv = Ainv
    def del_Ainv(self):
        """Delete the invert of A"""
        if hasattr(self, '_Ainv'):
            del self._Ainv

    Ainv = property(get_Ainv, set_Ainv, del_Ainv, doc='Invert of A')


    def interp(self, xo, yo, geterr=False, blockr=None):
        """Interpolate to positions xo,yo

        :Params:

            - **xo/yo**: Output positions.
            - **geterr**, optional: Also return errors.

        :Return: ``zo`` or ``zo,eo``
        """

        # Inits
        xo = N.asarray(xo, 'd')
        yo = N.asarray(yo, 'd')
        npo = xo.shape[0]
        vgf = self.variogram_func
        so = (self.nt, npo) if self.nt else npo
        zo = N.zeros(so, 'd')
        if geterr:
            eo = N.zeros(npo, 'd')
        if self.ncloud>1 or geterr:
            wo = N.zeros(npo, 'd')

        # Loop on clouds
        Ainv = self.Ainv
        next = int(not self._simple)
        for ic in xrange(self.ncloud): # TODO: multiproc here?

            # Distances to output points
            # dd = cdist(N.transpose([xi,yi]),N.transpose([xo,yo])) # TODO: test cdist
            dd = get_distances(xo, yo, self.xc[ic], self.yc[ic], mode=self.distfunc)

            # Form B
            np = self.npc[ic]
            B = N.empty((np+next, npo))
            B[:self.npc[ic]] = vgf(dd)
            if not self._simple:
                B[-1] = 1
            if self.exact:
                B[:np][isclose(B[:np], 0.)] = 0.
            del dd

            # Block kriging
            if blockr:
                tree = cKDTree(N.transpose([xo, yo]))
                Bb = B.copy()
                for i, iineigh in enumerate(tree.query_ball_tree(tree, blockr)):
                    Bb[:, i] = B[:, iineigh].mean()
                B = Bb

            # Compute weights
            W = N.ascontiguousarray(symm(Ainv[ic], N.asfortranarray(B, 'd')))

            # Simple kriging with adjusted mean for long distance values
            if self._simple and self.farvalue is not None:
                s = self.get_sill()
                Ais = self.get_sill() * Ainv[ic].sum(axis=0)
                mean = self.farvalue - (self.zc[ic] * Ais).sum()
                mean /= (1 - Ais.sum())
            else:
                mean = self.mean

            # Interpolate
            z = N.ascontiguousarray(dgemv(N.asfortranarray(W[:np].T, 'd'),
                N.asfortranarray(self.zc[ic]-mean, 'd')))
            if self._simple:
                z += mean

            # Simplest case
            if not geterr and self.ncloud<2:
                zo[:] = z.T
                continue

            # Get error
#            e = (W[:-1]*B[:-1]).sum(axis=0)
            e = (W*B).sum(axis=0)
            del W, B

            # Weigthed contribution based on errors
            w = 1/e**2
            if self.ncloud>1:
                z[:] *= w
            wo += w
            del w
            zo[:] += z.T ; del z

        # Error
        if geterr:
            eo = 1/N.sqrt(wo)

        # Normalization
        if self.ncloud>1:
            zo[:] /= wo

        gc.collect()
        if geterr: return zo, eo
        return zo

    __call__ = interp


class OrdinaryCloudKriger(CloudKriger):
    """Ordinary kriger using cloud splitting"""
    def __init__(self, x, y, z, mtype=None, vgf=None, npmax=1000,
            nproc=None, exact=False, distfunc='simple', errfunc=None,
            **kwargs):
        CloudKriger.__init__(self, x, y, z, 'ordinary', mtype=mtype, vgf=vgf, npmax=npmax,
            nproc=nproc, exact=False, distfunc=distfunc, errfunc=errfunc,
            **kwargs)
OrdinaryKriger = OrdinaryCloudKriger


class SimpleCloudKriger(CloudKriger):
    """Simple kriger using cloud splitting"""
    def __init__(self, x, y, z, mtype=None, vgf=None, npmax=1000,
            nproc=None, exact=False, distfunc='simple', errfunc=None,
            mean=None, farvalue=None, **kwargs):
        CloudKriger.__init__(self, x, y, z, 'simple', mtype=mtype, vgf=vgf, npmax=npmax,
            nproc=nproc, exact=False, distfunc=distfunc, errfunc=errfunc,
            mean=mean, farvalue=farvalue, **kwargs)


def krig(xi, yi, zi, xo, yo, vgf=None, geterr=False, **kwargs):
    """Quickly krig data"""
    return OrdinaryKriger(xi, yi, zi, vgf=vgf, **kwargs)(xo, yo, geterr=geterr)

def gauss3(x, y,
    x0=-1, y0=0.5, dx0=1, dy0=1, f0=1.,
    x1=1, y1=1, dx1=2, dy1=0.5, f1=-1,
    x2=0, y2=-1.5, dx2=.5, dy2=.5, f2=-.3,
    **kwargs):
    """Create data sample as function position and 3-gaussian function"""
    g = bivariate_normal(x, y, dx0, dy0, x0, y0)*f0
    g+= bivariate_normal(x, y, dx1, dy1, x1, y1)*f1
    g+= bivariate_normal(x, y, dx2, dy2, x2, y2)*f2
    g *= 10.
    return g

def gridded_gauss3(nx=100, ny=100, xmin=-3, xmax=3, ymin=-3, ymax=3, mesh=False, **kwargs):
    """Create a data sample on a grid using :func:`gauss3`"""
    x = N.linspace(xmin, xmax, nx)
    y = N.linspace(ymin, ymax, ny)
    xx, yy = N.meshgrid(x, y)
    zz = gauss3(xx, yy, **kwargs)
    if mesh:
        return xx, yy, zz
    return x, y, zz

def random_gauss3(**kwargs):
    """Create a data sample of random points using :func:`gauss3`"""
    x, y = random_points(**kwargs)
    z = gauss3(x, y, **kwargs)
    return x, y, z

def random_points(np=200, xmin=-3, xmax=3, ymin=-3, ymax=3, **kwargs):
    """Generate random coordinates of points"""
    x = P.rand(np)*(xmax-xmin)+xmin
    y = P.rand(np)*(ymax-ymin)+ymin
    return x, y




