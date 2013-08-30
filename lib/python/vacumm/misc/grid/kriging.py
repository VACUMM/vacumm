# -*- coding: utf8 -*-
"""
Kriging utilities inspired from the AMBHAS library (http://www.ambhas.com/).
"""
# Copyright or © or Copr. Actimar (contributor(s) : Stephane Raynaud) (2013)
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

import multiprocessing
from multiprocessing import Pool,  cpu_count

import numpy as N
import pylab as P
from scipy.optimize import curve_fit
import scipy.linalg.fblas
blas_dgemv = scipy.linalg.fblas.dgemv
def dgemv(a, x): return blas_dgemv(1., a, x)
try:
    from _blaslapack import symm, sytri
except:
#    print 'Falling back to builtin functions'
    dgemm = scipy.linalg.fblas.dgemm
    def symm(a, b): return dgemm(1., a, b)
    sytri = N.linalg.inv
import gc

class KrigingError(Exception):
    pass

#: Variogram model types
variogram_model_types = ['linear', 'exponential', 'spherical', 'gaussian']

def variogram_model_type(mtype=None):
    """Check the the variogram model type
    
    :Params: 
    
        - **mtype**, optional: ``None``, and index or a string matching
          an element of :data:`variogram_model_types`.
          If set to ``None``, it defaults to ``"exponential"``.
    """
    if mtype is True:
        return variogram_model_types
    errmsg = []
    for i, vtype in enumerate(variogram_model_types):
        errmsg += '"%s" (=%i)'%(vtype, i)
    errmsg = 'Invalid variogram model type. Please choose one of: '+ ', '.join(errmsg)
    if mtype is None: mtype = 'exponential'
    if isinstance(mtype, int):
        if i<0 or i>len(variogram_model_types)-1:
            raise KrigingError(errmsg)
        return variogram_model_types[i]
    if not isinstance(mtype, basestring):
        raise KrigingError(errmsg)
    for vtype in variogram_model_types:
        if vtype.startswith(mtype): return vtype
    raise KrigingError(errmsg)
        

def variogram_model(mtype, n, s, r=None, nrelmax=0.2):
    """Get the variogram model function from its name"""
    
    mtype = variogram_model_type(mtype)
    
    n = max(n, 0)
    
    if mtype == 'linear': 
        return lambda h: n + h * s # FIXME: h/r?
        
    n = min(n, nrelmax*s)
    if r is None:
        raise KrigingError('Please specifiy the range parameter (r) for the variogram model')

    if mtype=='exponential':
        return lambda h: n + (s-n) * (1 - N.exp(-3*h/r))
        
    if mtype=='spherical':
        return lambda h: n + (s-n)*((1.5*h/r - 0.5*(h/r)**3)*int(h<=r) + int(h>r))
        
    if mtype=='gaussian':
        return lambda h: n + (s-n)*(1-N.exp(-3*h**2/r**2))
        
    raise KrigingError(errmsg)

class VariogramModel(object):
    def __init__(self, mtype, fpar=None):
        self.mtype = variogram_model_type(mtype)
        self.fpar = fpar
    @staticmethod
    def fix_params(args, fpar):
        if hasattr(fpar, '__len__'):
            args = list(args)
            fpar = fpar[:len(args)+1]
            for i, farg in enumerate(fpar):
                if farg is not None:
                    args[i] = farg     
        return args   
    def __call__(self, d, *args):
        args = self.fix_params(args, self.fpar)
        return variogram_model(self.mtype, *args)(d)

def _get_xyz_(x, y, z, check=True):
    if not x.ndim==1 and x.shape!=y.shape and x.shape!=z.shape:
        raise KrigingError('x, y and z must have the same 1D shape')
    mask = N.ma.nomask
    if N.ma.isMA(z): mask = N.ma.getmaskarray(z)   
    if N.ma.isMA(x): mask |= N.ma.getmaskarray(x)
    if N.ma.isMA(y): mask |= N.ma.getmaskarray(y)
    if mask is not N.ma.nomask and mask.any():
        good = ~mask
        if check and not good.any():
            raise KrigingError('All data are masked')
        try:
            x = x.compress(good)
        except:
            pass
        y = y.compress(good)
        z = z.compress(good)
    return x, y, z
    
def variogram(x, y, z, binned=None, nmax=1500, nbindef=30, nbin0=None, nbmin=10):
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
          
    """
    x, y, z = _get_xyz_(x, y, z)
    npts = x.shape[0] 
    
    # Undepsample?
    if npts>2*nmax:
        samp = npts/nmax
        x = x[::samp]
        y = y[::samp]
        z = z[::samp]
     
    # Distances
    x0, x1 = N.meshgrid(x, x)
    y0, y1 = N.meshgrid(y, y)
    dd = N.sqrt((x1-x0)**2+(y1-y0)**2)
    del x0, x1, y0, y1
    
    # Variogram
    z0, z1 = N.meshgrid(z, z)
    vv = 0.5*(z1-z0)**2
    del z0, z1

    # Unique
    iup = N.triu_indices(dd.shape[0])
    d = dd[iup]
    v = vv[iup]
    del dd, vv
    
    # Direct variogram?
    if binned is None and len(d)>nbindef*nbmin: binned = True
    if binned is True: binned = nbindef
    if not binned: return d, v
    
    # Rebin
    ii = N.argsort(d)
    np = d.shape[0]
    nbin = binned
    bins = N.linspace(0, np-1, binned+1).astype('i').tolist()
    if nbin0 is None:
        if bins[1]>2*nbmin:
            nbin0 = min(bins[1]/nbmin, binned)
        else:
            nbin0 = 0
    if nbin0>1: # split first bin
        bins = N.linspace(0., bins[1], nbin0+1).astype('i')[:-1].tolist()+bins[1:]
        nbin = len(bins)
    db = N.empty(nbin)
    vb = N.empty(nbin)
    for ib in xrange(nbin-1):
        iib = ii[bins[ib]:bins[ib+1]+1]
        db[ib] = d[iib].mean()
        vb[ib] = v[iib].mean()
    return db, vb
    
def variogram_fit(x, y, z, mtype, getp=False, geterr=False, fpar=None):
    """Fit a variogram model to data and return the function"""
    # Estimated variogram
    d, v = variogram(x, y, z)
    
    # Variogram model
    vm = VariogramModel(mtype, fpar=fpar)
    
    # First guess of paramaters
    imax = N.argmax(v)
    p0 = [0., v[imax]]
    if mtype!='linear': 
        p0.append(d[imax])
    
    # Fitting
    p, e = curve_fit(vm, d, v, p0=p0)
    res = p if getp else variogram_model(vm.mtype, *p)
    if not geterr: return res
    return res,  (variogram_model(vm.mtype, *p)(d)-v).std()
    
def variogram_multifit(xx, yy, zz, mtype=None, getp=False, fpar=None):
    """Same as :func:`variogram_fit` but with several samples"""
    pp = []
    for i, (x, y, z) in enumerate(zip(xx, yy, zz)):
        x, y, z = _get_xyz_(x, y, z, check=False)
        if len(x)==0: continue
        p = variogram_fit(x, y, z, mtype, getp=True)
        if fpar is not None: 
            p = VariogramModel.fix_params(p, fpar)
        pp.append(p)
    pp = N.asarray(pp)
    if pp.shape[0]==0:
        raise KrigingError('All data are masked')
    mp = N.median(pp, axis=0)
    return variogram_model(mtype, *mp)

def cloud_split(x, y, npmax=1000, getdist=True, getcent=True):
    """Split data intot cloud of points of max size cmax:
    
    :Returns: ``None`` if ``len(x)<=cmax``
    
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
    res = sytri(N.asfortranarray(A, 'd'))
    if isinstance(res, tuple):
        info = res[1]
        if info: raise KrigingError('Error during call to Lapack DSYTRI (info=%i)'%info)
        return res[0]
    else:
        return res

class OrdinaryKriger(object):
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
           Weights are estimated using the max
           of ``1/B`` where ``B`` is output
           variogram matrix.
    
    :Params:
    
        - **x/y/z**: Input positions and data (masked array).
        - **mtype**, optional: Variogram model type (defaults to 'exp').
          See :func:`variogram_model_type` and :func:`variogram_model_types`.
        - **vgf**, optional: Variogram function. If not set,
          it is estimated using :meth:`variogram_fit`.
        - **npmax**, optional: Maxima size of cloud.
        - **nproc**, optional: Number of processes to use
          to invert matrices. Set it to a number <2 to switch off
          parallelisation.
          
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
    
    def __init__(self, x, y, z, mtype=None, vgf=None, npmax=1000, nproc=None):
        self.x, self.y, self.z = _get_xyz_(x, y, z)
        self.np = self.x.shape[0]
        self.mtype = variogram_model_type(mtype)
        self.npmax = npmax
        if nproc is None: 
            self.nproc = cpu_count()
        else:
            self.nproc = max(1, min(cpu_count(), nproc))
        self._setup_clouds_()
        if callable(vgf): 
            self.variogram_func = vgf
            
    def __len__(self):
        return self.x.shape[0]
    
    def _get_xyz_(self, x=None, y=None, z=None):
        if x is None: x = self.x
        if y is None: y = self.y
        if z is None: z = self.z
        return _get_xyz_(x, y, z)
        
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
                    getattr(self, xyz+'c').append(getattr(self, xyz)[indices[ic]])
                
                # Size
                self.npc.append(len(indices[ic]))
                
        else: # Single cloud
            
            self.xc = [self.x]
            self.yc = [self.y]
            self.zc = [self.z]
            self.ncloud = 1
            self.npc = [self.np]
            self.cwfunc = [lambda x, y: 1.]
            
    def plot_clouds(self, **kwargs):
        """Quickly Plot inputs points splitted in clouds"""
        P.figure()
        for x, y in zip(self.xc, self.yc):
            P.plot(x, y, **kwargs)
        P.show()
            
    def variogram_fit(self, x=None, y=None, z=None):
        """Estimate the variogram function by using :func:`variogram_fit`"""
        x, y, z = self._get_xyz_(x, y, z)
        self.variogram_func = variogram_fit(x, y, z, self.mtype)
        return self.variogram_func
        
    def variogram_multifit(self, xx, yy, zz):
        """Estimate the variogram function by using :func:`variogram_multifit`"""
        del self.variogram_func
        self.variogram_func = variogram_multifit(xx, yy, zz, self.mtype)
        return self.variogram_func
        
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
        """Get the invert of A"""

        # Variogram function
        vgf = self.variogram_func
        
        # Already computed
        if hasattr(self, '_Ainv'): 
            return self._Ainv
        
        # Loop on clouds
        if not hasattr(self, '_dd'): self._dd = []
        Ainv = []
        AA = []
        for ic in xrange(self.ncloud):
            
            # Get distance between input points
            if len(self._dd)<ic+1:
                x0, x1 = N.meshgrid(self.xc[ic], self.xc[ic])
                y0, y1 = N.meshgrid(self.yc[ic], self.yc[ic])
                dd = (x1-x0)**2 ; del x0,x1
                dd += (y1-y0)**2 ; del y0,y1
                dd = N.sqrt(dd).astype('d')
                self._dd.append(dd)
            else:
                dd = self._dd[ic]
                
            # Form A
            A = N.empty((self.npc[ic]+1, self.npc[ic]+1))
            A[:-1, :-1] = vgf(dd)
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
            Ainv = Pool(self.nproc).map(syminv, AA, chunksize=1)
        
        # Fortran arrays
        Ainv = [N.asfortranarray(ainv, 'd')for ainv in Ainv]
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


    def interp(self, xo, yo, geterr=False):
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
        zo = N.zeros(npo, 'd')
        if geterr: eo = zo.copy()
        if self.ncloud>1:
            wo = zo.copy()
        
        # Loop on clouds
        Ainv = self.Ainv
        for ic in xrange(self.ncloud): # TODO: multiproc here?
        
            # Distances to output points
            xxo, xxi = N.meshgrid(xo, self.xc[ic])
            dd = (xxo-xxi)**2 ; del xxi, xxo
            yyo, yyi = N.meshgrid(yo, self.yc[ic])
            dd += (yyo-yyi)**2 ; del yyi, yyo
            dd = N.sqrt(dd)
            
            # Form B
            B = N.empty((self.npc[ic]+1, npo))
            B[-1] = 1
            B[:-1] = vgf(dd)
            del dd
            
            # Compute weights
            W = N.ascontiguousarray(symm(Ainv[ic], N.asfortranarray(B, 'd')))
            
            # Interpolate
            z = N.ascontiguousarray(dgemv(N.asfortranarray(W[:-1].T, 'd'), 
                N.asfortranarray(self.zc[ic], 'd')))
                
            # Weigthed contribution
            if self.ncloud>1:
                w = 1/B[:-1].min(axis=0)
                z[:] *= w
                wo += w
                if not geterr: del w
            zo[:] += z ; del z
                        
            # Get error
            if geterr: 
                e = (W[:-1]*B[:-1]).sum(axis=0)
                if self.ncloud>1: 
                    e*= w
                    del w
                eo[:] += e
            
            del W, B
            
        # Normalization
        if self.ncloud>1:
            zo[:] /= wo
            if geterr:
                eo[:] /= wo
            del wo
        
        gc.collect()
        if geterr: return zo, eo
        return zo
    
    __call__ = interp

def krig(xi, yi, zi, xo, yo, vgf=None, geterr=False, **kwargs):
    """Quickly krig data"""
    return OrdinaryKriger(xi, yi, zi, vgf=vgf, **kwargs)(xo, yo, geterr=geterr)

def gauss3(x, y, 
    x0=-1, y0=0.5, dx0=1, dy0=1, f0=1., 
    x1=1, y1=1, dx1=2, dy1=0.5, f1=-1, 
    x2=0, y2=-1.5, dx2=.5, dy2=.5, f2=-.3, 
    **kwargs):
    """Create data sample as function position and 3-gaussian function"""
    g = P.bivariate_normal(x, y, dx0, dy0, x0, y0)*f0
    g+= P.bivariate_normal(x, y, dx1, dy1, x1, y1)*f1
    g+= P.bivariate_normal(x, y, dx2, dy2, x2, y2)*f2
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
    

    
    
