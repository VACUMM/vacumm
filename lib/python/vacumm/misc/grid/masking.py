# -*- coding: utf8 -*-
"""Utilities to deal with masks and selections

.. seealso::

    Tutorials: :ref:`user.tut.misc.grid.masking`, :ref:`user.tut.misc.grid.polygons`
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
#print 'importing masking 0'

import numpy as N
MA = N.ma
from cdms2.hgrid import TransientCurveGrid
from cdms2.grid import AbstractGrid
import cdms2, MV2
from mpl_toolkits.basemap import Basemap
from _geoslib import Point, Polygon, LineString
import gc
from genutil import minmax
cdms = cdms2

__all__ = ['get_coast', 'get_coastal_indices', 'GetLakes', 'polygon_mask', 'masked_polygon', 'polygons', 'd2m', 'polygon_select', 'envelop', 'check_poly_islands', 'check_poly_straits','t2uvmasks', 'mask2d', 'grid_envelop', 'convex_hull', 'uniq', 'rsamp', 'zcompress', 'Lakes', 'erode_coast', 'resol_mask', 'get_dist_to_coast']
__all__.sort()

def get_coast(mask, land=True, b=True, borders=True, corners=True):
    """Get a mask integer of ocean coastal points from a 2D mask

    - **mask**: Input mask with 0 on water and 1 on land.
    - *land*: If True, coast is on land [default: True]
    - *corners*: If True, consider borders as coastal points [default: True]
    - *borders:* If True, consider corners as coastal points [default: True]
    
    At each point, return 0 if not coast, else an interger ranging from 1 to (2**8-1) to describe the coast point::
    
        128 4  64
        8   +   2
        16  1  32
    """
    # Initialize mask and its boundaries
    sh = list(mask.shape)
    imask = N.zeros((sh[0]+2, sh[1]+2), 'i')
    imask[1:-1, 1:-1] = mask.astype('i')
    imask[0] = imask[1]
    imask[-1] = imask[-2]
    imask[:, 0] = imask[:, 1]
    imask[:, -1] = imask[:, -2]
    for j in 0, -1:
        jp = N.sign(j)
        for i in 0, -1:
            ip = N.sign(i)
            imask[j, i] = imask[j+jp, j] + imask[j, i+ip]
    # Coast on land or ocean
    if not land:
        iimask = imask
        imask = 1-imask
    else:
        iimask = 1-imask
    # Coast along borders only
    mborders = \
        1  *iimask[ :-2,1:-1]+ \
        2  *iimask[1:-1,2:  ]+ \
        4  *iimask[2:  ,1:-1]+ \
        8  *iimask[1:-1, :-2]
    # Coast as corners only
    mcorners = \
        16 *iimask[ :-2, :-2]+ \
        32 *iimask[ :-2,2:  ]+ \
        64 *iimask[2:  ,2:  ]+ \
        128*iimask[2:  , :-2]

    mcorners = N.where(N.equal(mborders, 0), mcorners, 0)
    m = imask[1:-1,1:-1]*(borders*mborders+corners*mcorners)
    del imask, iimask
    gc.collect()
    # Return
    if b: return m.astype(bool)
    return m

def get_coastal_indices(mask, coast=None, asiijj=False, **kwargs):
    """Get indices of ocean coastal points from a 2D mask

    :Params:
    
        - **mask**: Input mask with 0 on water and 1 on land.
        - *boundary*: If True, consider outside boundary as land [default: True]
        - *asiijj*: Get results as II,JJ instead of [(j0,i0),(j1,i1)...]
    """
    if mask is MA.nomask: return []
    if coast is None:
        coast = get_coast(mask, **kwargs)
    iy,ix = N.indices(mask.shape)
    jj = N.compress(coast.ravel(),iy.ravel())
    ii = N.compress(coast.ravel(),ix.ravel())
    if asiijj: return ii, jj
    return zip(jj, ii)
    
def get_dist_to_coast(grid, mask=None, proj=True):
    """Get the distance to coast on grid
    
    :Params:
    
        - **grid**: (x,y), grid or variable with a grid.
        - **mask**: Land/sea mask or variable with a mask. If not provided, 
          gets mask from ``grid`` if it is a masked variable,
          or estimates it using :func:`polygon_mask`
          and a GHSSC shoreline resolution given by
          :func:`~vacumm.misc.grid.basemap.gshhs_autores`.
        - **proj**: Distance are computed by converting coordinates
          from degrees to meters.
          
    :Example:
    
        >>> dist = get_dist_to_coast(sst.getGrid(), sst.mask)
        >>> dist = get_dist_to_coast(sst)
        
    :See also: :func:`get_coastal_indices`
    """
    if mask is None and N.ma.isMA(grid):
        mask = grid
    xx, yy = M.get_xy(grid, proj=proj, mesh=True, num=True)
    if mask is None:
        from basemap import gshhs_autores, merc
        x, y = M.get_xy(grid, proj=False, num=True, mesh=False)
        res = gshhs_autores(x.min(), x.max(), y.min(), y.max())
        mask = polygon_mask(grid, res)
    elif N.ma.isMA(mask):
        mask = N.ma.getmaskarray(mask)
    mask = mask2d(mask)
    ii, jj = get_coastal_indices(mask, asiijj=True)
    good = ~mask 
    xc = xx[jj, ii]
    yc = yy[jj, ii]
    nc = len(xc)
    ng = good.sum()
    xdiff = N.resize(xx[good], (nc, ng)) - N.resize(xc, (ng, nc)).T
    ydiff = N.resize(yy[good], (nc, ng)) - N.resize(yc, (ng, nc)).T
    dists = xdiff**2+ydiff**2 ; del xdiff, ydiff
    md = N.sqrt(dists.min(axis=0)) ; del dists
    mindists = N.ma.zeros(xx.shape)
    mindists[mask] = N.ma.masked
    mindists[good] = md
    del good
    if cdms2.isVariable(grid):
        mindists = MV2.asarray(mindists)
        mindists.id = 'coastdist'
        mindists.long_name = 'Distance to coast'
        mindists.units = 'm'
        M.set_grid(mindists, M.get_grid(grid))
    return mindists
    
    

def erode_coast(var, mask2d=None, copy=True, maxiter=10, corners=0., hardmask=None, **kwargs):
    """Erode coastal mask of ``var`` to fit to ``mask2d`` using 2D smoothing
    
    Soothing is performed using :func:`~vacumm.misc.filters.shapiro2d`, 
    in a convergence loop. Loop is ended if :
    
        - the mask no longer changes,
        - the number of iteration exceeds ``maxiter``.
    
    .. warning:: 
    
        Must be used only for erodeing a thin border
        of the coast.
    
    :Params:
    
        - **var**: Atleast-2D A :mod:`MV2` variable with mask.
        - **mask2d**, optional: Reference 2D mask.
        - **maxiter**, optional: Maximal number of iteration for convergence loop. 
          If equal to ``None``, no check is performed.
        - **corners**, optional: Weights of corners when calling :func:`~vacumm.misc.filters.shapiro2d`.
        - **copy**, optional: Modify the variable in place or copy it.
        - **hardmask**, optional: Mask always applied after an erosion step.
        - Other keywords are passed to :func:`~vacumm.misc.filters.shapiro2d`.
        
    :Return: 
    
        - A :mod:`MV2` variable
        
    :Tutorial: :ref:`user.tut.misc.grid.masking.erode_coast`
    """
    from vacumm.misc.filters import shapiro2d
    varo = var.clone() if copy else var
    if mask2d is None:
        mask2d = N.zeros(var.shape[-2:], '?')
    masknd = mask2d if (var.ndim==2) else N.resize(mask2d, var.shape)
    if hardmask is not None and hardmask.ndim != var.ndim:
        hardmask = N.resize(hardmask, var.shape)
    # Not masked
    if varo.mask is MV2.nomask or not var.mask.any():
        varo[:] = MV2.masked_where(masknd, varo)
        return varo
    # Convergence loop
    oldnmasked = varo.mask.sum()
    i=0
    while maxiter is None or i<maxiter:
        varc = shapiro2d(varo, corners=corners, copy=True, mask='minimal', **kwargs)
        varo[:] = MV2.where(varo.mask, varc, varo)
        if hardmask is not None:
            varo[:] = MV2.masked_where(hardmask, varo, copy=0)
        del varc
        maskc = varo.mask | masknd
        nmasked = maskc.sum()
        if nmasked == oldnmasked: break
        del maskc
        oldnmasked = nmasked
        i+=1
    varo[:] = MV2.masked_where(masknd, varo, copy=0)
    return varo
        

class Lakes(object):
    """Find lakes in a 2D mask where 0 is water and 1 is land

    :Example:
    
        >>> from vacumm.misc.grid import Lakes
        >>> import numpy as N
        >>> from pylab import pcolor,show,title
        >>> # Build the mask
        >>> mask=N.ones((500,500))
        >>> mask[3:10,100:102]=0
        >>> mask[103:110,200:210]=0
        >>> mask[203:210,200:220]=0
        >>> # Find the lakes
        >>> lakes = Lakes(mask,nmaxcells=80)
        >>> print 'Number of lakes:', len(lakes.indices())
        Number of lakes: 2
        >>> pcolor(lakes.mask(),shading=False)
        >>> title('Lakes')
        >>> show()
        >>> first_two_lakes = lakes[:2]
    """
    def __init__(self,mask,nmaxcells=None):

        # Inits
        self._nmaxcells = nmaxcells
        self._ny,self._nx = mask.shape
        self._mask = N.array(mask[:], dtype='?')
        self._imask = self._mask.astype('i')
        self._iy, self._ix = N.indices(mask.shape)
        self._np = self._imask.sum()
        
        # Find lakes using labelling
        import scipy.ndimage
        self._lakes, self.nlakes = scipy.ndimage.label(1-self._imask)
        self._ncells = [N.equal(self._lakes, ilake).astype('i').sum() for ilake in xrange(1, self.nlakes+1)]
    
        # Compute indices
        self._indices = []
        for ilake in xrange(1, self.nlakes+1):
            mask = self._lakes == ilake
            self._indices.append(zip(self._ix[mask], self._iy[mask]))
        
        # Sorting argument
        self._argsort = N.argsort(self._ncells)[::-1]


    def __call__(self,**kwargs):
        return self.indices(**kwargs)
        
    def __len__(self):
        return self.nlakes
        
    def __getitem__(self, item):
        lakes = N.zeros(self._lakes.shape, 'l')
        for lakeid in self._size2ids_(item):
            lakes[self._lakes==lakeid] = lakeid
        return lakes
            
    def _size2ids_(self, item):
        ids = N.arange(self.nlakes)[self._argsort[item]]
        ids += 1
        if isinstance(ids, int): ids = N.asarray([ids])
        return ids

    def plot(self, **kwargs):
        """Display the lakes
        
        Keywords are passed to :func`~matplotlib.pyplot.pcolor`
        """
        ny, nx = self._mask.shape
        xx, yy = N.meshgrid(N.arange(nx+1), N.arange(ny+1))
        import pylab as P
        P.pcolor(xx, yy, self.lakes(), **kwargs)
        P.show()
    
    def ncells(self):
        """Number of cell point for each lake"""
        return self._ncells


    def indices(self, id=None, nbig=None):
        """Return a list of lake coordinates
        
        - *id*: Select lake #<id>
        - *nbig*: Consider the first ``nbig`` greatest lakes
        """
        if id is not None:
            return self._indices[id]
        if nbig is None:
            return self._indices
        indices = []
        for lakeid in self._size2ids_(slice(0, nbig)):
            indices.append(self._indices[lakeid])
        return indices

    def lakes(self, id=None, nbig=None):
        """Return an array of same size as masks, but with points in lakes set to its lake id, and others to set to zero
        
        - *id*: Select lake #<id>
        - *nbig*: Consider the first ``nbig`` greatest lakes
        """
        if id is None:
            if nbig is not None:
                return self[:nbig]
            return self._lakes.copy()
        return N.where(self._lakes==id, self._lakes, 0)

    def mask(self, id=None, nbig=None, **kwargs):
        """Returns a boolean land/sea mask from lakes() (land is True)
        
        - *id*: Select lake #<id>
        - *nmaxcells*: Do not consider lakes with number of points greater than nmaxcells
        
        :Example:
        
        >>> mask_lake2 = GetLakes(mask).mask(2)
        """
        return self.lakes(id=id, nbig=nbig) == 0
        
    def ocean(self, mask=False):
        """Get the biggest lake or its mask (integer type)
        
        - *mask*: If True, return the mask, not the lake [default: True]
        """
        if mask:
            return self.mask(nbig=1)
        return self[0]

        
GetLakes = Lakes

def polygon_mask(gg, polys, mode='intersect', thresholds=[.5, .75], ocean=False, fractions=0, yclean=True, premask=None):
    """Create a mask on a regular or curvilinear grid according to a polygon list
    
    :Params:
    
        - **gg**: A cdms grid or variable or a tuple of (xx,yy).
        - **polys**: A list of polygons or GMT resolutions or Shapes instance like shorelines.
        - **mode**, optional: Way to decide if a grid point is masked. Possible values are:
        
            - ``intersect``, 1, ``area`` (default): Masked if land area fraction is > ``thresholds[0]``. If more than one intersections, leand area fraction must be > ``thresholds[1]`` to prevent masking straits.
            - else: Masked if grid point inside polygon.
            
        - **thresholds**, optional: See ``intersect`` mode [default: [.5, .75]]
        - **ocean**, optional: If True, keep only the ocean (= biggest lake) [default: True]
        
    :Result:
    
        A :mod:`numpy` array of boolean (mask) with ``False`` on land.
    """
    
    # Get proper numeric axes
    gg = M.get_grid(gg)
    curved = M.isgrid(gg, curv=True)
    xx, yy = M.get_xy(gg, mesh=True, num=True)
    
    # Thresholds
    if not isinstance(thresholds, (list, tuple)):
        thresholds = [thresholds, thresholds]
    elif len(thresholds) == 1:
        thresholds += thresholds
            
    # Bounds if mode is 'intersect'
    mask = N.zeros(xx.shape, dtype='?')
    bmode = mode in ['intersect', 1, 'area']
    fmode = mode in [2, 'fractions', 3]
    if fmode: bmode = True
    xxb = M.bounds2d(xx)
    yyb = M.bounds2d(yy)
    if bmode:
#        yyb[yyb>90.] = 89.99
#        yyb[yyb<-90.] = -89.99
        dx = xxb.ptp()/xxb.shape[1]
        dy = yyb.ptp()/yyb.shape[0]
        clip = (xxb.min()-dx, yyb.min()-dy, xxb.max()+dx, yyb.max()+dy)
        ymins = yyb[:, 0, 0]-dy/10.
        xmin = xxb[0, 0].min()-dx/10.
        xmax = xxb[0, -1].max()+dx/10.
        ymax = yyb[-1, 0].max()+dy/10.
    else:
        dx = xx.ptp()/xx.shape[1]
        dy = yy.ptp()/yy.shape[0]
        clip = (xx.min()-dx, yy.min()-dy, xx.max()+dx, yy.max()+dy)
        ymins = yy[:, 0]-dy/10.
        xmin = xx[0, 0]-dx/10.
        xmax = xx[0, -1]+dx/10.
        ymax = yy[-1, 0]+dy/10.
    if curved:
        xxc, yyc = M.meshcells(xx, yy)
        dxx = N.diff(xxc, axis=1)
        dxy = N.diff(xxc, axis=0)
        dyx = N.diff(yyc, axis=1)
        dyy = N.diff(yyc, axis=0)
        for var in dxx, dxy, dyx, dyy:
            var /= 10.
        
    if premask is not None:
        premask = N.asarray(premask, 'i')
    
    # Get instances of Polygons
    polys = polygons(polys, clip=clip, shapetype=2)

    # Loop on grid points
    skipped = 0
    for j in xrange(xx.shape[0]):
        
        # Get polygons of the curent row
        ypolys = []
        if curved:
            # base
            rowdown = N.vstack((xxc[j], yyc[j])).T
            rowup = N.vstack((xxc[j+1], yyc[j+1])).T
            # bottom margin
            rowdown[1:-1, 0] -= .5*(dyx[j, :-1]+dyx[j, 1:])
            rowdown[1:-1, 1] -= dyy[j, 1:-1]
            # top margin
            rowup[1:-1, 0] += .5*(dyx[j, :-1]+dyx[j, 1:])
            rowup[1:-1, 1] += dyy[j, 1:-1]
            # margin at corners
            for ij, ddx, ddy in (0, dxx, dxy), (1, dyx, dyy): # x,y
                rowdown[0, ij] -= ddx[j, 0] + ddy[j, 0]
                rowdown[-1, ij] += ddx[j, -1] - ddy[j, -1]
                rowup[0, ij] -= ddx[j+1, 0] - ddy[j, 0]
                rowup[-1, ij] += ddx[j+1, -1] + ddy[j, -1]
            rowpoly_array = N.vstack((rowdown, rowup[::-1]))
            del rowdown, rowup
        else:
            rowpoly_array = N.array([
                [xmin, yyb[j, 0, 0]-dy/10.], [xmax, yyb[j, 0, 0]-dy/10.], 
                [xmax, yyb[j, 0, 3]+dy/10.], [xmin, yyb[j, 0, 3]+dy/10.]])
        rowpoly = Polygon(rowpoly_array)
        for poly in polys:
            try:
                if poly.intersects(rowpoly):
                    ypolys.extend(poly.intersection(rowpoly))
            except:
                pass
        del rowpoly_array, rowpoly

        for i in xrange(xx.shape[1]):

            # Pre-masking
            if premask is not None and premask[j, i] != -1:
                mask[j, i] = premask[j, i]
                skipped += 1
                continue
        
            # Point is inside
            if not bmode:
                for poly in ypolys:
                    mask[j, i] = Point((xx[j, i],yy[j, i])).within(poly)
                    if mask[j, i]: break
                continue
            
            # Check the cell
            cell_poly = Polygon(N.asarray([xxb[j, i, :], yyb[j, i, :]]).transpose())
            if isinstance(cell_poly, (LineString, Point)):
                raise TypeError, 'cell'
            intersect_area = 0.
            nintersect = 0
            for ip, poly in enumerate(ypolys):
                if poly is None: continue
                if isinstance(poly, (LineString, Point)):
                    raise TypeError, 'poly'
                    
                # Check X range
                if xxb[j, i, :].max() <= poly.boundary[:, 0].min()  or \
                    xxb[j, i, :].min() >= poly.boundary[:, 0].max(): 
                        continue
                    
                # Check if polygon covers cell
                if cell_poly.within(poly):
                    nintersect = 1
                    intersect_area = cell_poly.area()
                    break
                    
                # One polygon
                ok = False
                if cell_poly.intersects(poly):
                    
                    # Remove from list of polys if inside cell
                    if poly.within(cell_poly):
                        intersect_area += poly.area()
                        nintersect += 1
                        del ypolys[ip]
                        continue
                    try:
                        intersections = cell_poly.intersection(poly)
                    except:
#                        print '[polygon_mask] bad intersection polys', i, j
                        continue
                    for intersection in intersections:
                        if isinstance(intersection, (LineString, Point)):
                            continue
                        # One intersection polygon
                        this_area = intersection.area()
                        intersect_area += this_area
                        nintersect += 1
                        if this_area == poly.area():
                            ok = True
                            break
                    del intersections
                if ok: break
                
            # Fraction
            land_fraction = intersect_area/cell_poly.area()
            if fmode: # Fraction only
                mask[j, i] = land_fraction
                continue
                
            # Mask according to fraction (basic masking + strait checking)
            mask[j, i] = land_fraction >= thresholds[nintersect > 1]
#           mask[j, i] = intersect_fraction < thresholds[1] and nintersect > 1
            del cell_poly
            
        for ypoly in ypolys: del ypoly
        del ypolys
#   print 'skipped', skipped
    # Select ocean only
    if ocean:
        return GetLakes(mask).ocean()
    
    return mask




def masked_polygon(vv, polys, copy=0, **kwargs):
    """Mask a variable according to a polygon list
    
    - **vv**: The variable
    - **polys**: Polygons (or something like that)
    - *copy*: Copy data [default: 0]
    - Other keywords are passed to :func:`MV.masked_where`
    """
#   # Get grid
#   gg = M.get_xy(vv, m=None)
#   
    # Get mask
    mask = polygon_mask(vv, polys, **kwargs)
    
    # Mask the variable
    vv2 = vv.clone() if copy else vv
    vv2[:] = MV2.masked_where(mask, vv, copy=0)
    del mask
    return vv2

def masked_ocean(vv, polys=None, **kwargs):
    """Remove lakes from a variable to keep only ocean where ocean the the biggest lake
    
    - **vv**: The variable
    """
    pass


def polygons(polys, proj=None, clip=None, shapetype=2, **kwargs):
    """Return a list of Polygon instances
    
    - **polys**: A tuple or list of polygon proxies (see examples).
    - *shapetype*: 1 = Polygons, 2=Polylines(=LineStrings) [default: 2]
    - *proj*: Geographical projection to convert positions.
    
    Example :
    
    >>> from vacumm.misc.grid.masking import polygons
    >>> import numpy as N
    >>> X = [0,1,1,0]
    >>> Y = N.array([0,0,1,1])
    >>> polygons( ([X,Y],) )
    >>> polygons( (zip(X,Y), [X,Y], N.array([X,Y]) )
    >>> polygons( (polygons(([X,Y],), polygons(([X,Y],)))
    >>> polygons( ([min(X),min(Y),max(X),max(Y)],) )
    """
    
    # Input
    if isinstance(polys, (Basemap, str, Polygon, N.ndarray, Shapes, AbstractGrid, tuple)):
        polys = [polys]
        
    if kwargs.has_key('m'): proj = kwargs['m']
        
    # Numeric or geos polygons
    out_polys = []
    gmt_polys = []
    if shapetype == 1:
        shaper = LineString
    else:
        shaper = Polygon
    if clip is not None: clip = polygons([clip], proj=proj)[0]
    from misc import isgrid, isrect, curv2rect
    for poly in polys:
        
        # Grid (argument to get_grid)
        if cdms2.isVariable(poly) or isgrid(poly) or \
            (isinstance(poly, tuple) and len(poly)==2):
            grid = M.get_grid(poly)
            if grid is None: continue
            poly = grid_envelop(grid)
            
        # It's already a polygon
        if isinstance(poly, Polygon): 
            if clip is not None:
                if poly.intersects(clip):
                    try:
                        out_polys.extend(poly.intersection(clip))
                    except:
                        pass
            else:
                out_polys.append(poly)
            continue
        
        # Polygon from GMT
        if isinstance(poly, (str, Basemap)):
#           gmt_polys.append(poly)
            out_polys.extend(GSHHS_BM(poly, clip=clip))
            continue
            
        # Shapes instance
        if isinstance(poly, Shapes):
            # Good polygon?
            if poly.get_type() == 0 or poly.get_type() != shapetype:
                continue
            # Clip
            poly.clip(clip)
            # Append to list
            out_polys.extend(poly.get_shapes())
            continue
        
        # Convert to array
        poly = N.asarray(poly, 'float64')
        
        # xmin,ymin,xmax,ymax form
        if poly.ndim == 1:
            assert len(poly) == 4, '1D form must have 4 elements (xmin,ymin,xmax,ymax), not %i'%len(poly)
            xmin,ymin,xmax,ymax = poly
            poly = N.asarray([[xmin, ymin], [xmax, ymin], [xmax, ymax], [xmin, ymax]])
            
        # Check order and geographic projection
        if poly.shape[0] == 2:
            if callable(proj):
                poly[:] = proj(*poly)
            poly = poly.transpose()
        elif callable(proj):
            poly[:, 0], poly[:, 1] = proj(poly[:, 0], poly[:, 1])
    
        # Get polygon object
        poly = shaper(poly)
        if clip is not None:
            if poly.intersects(clip):
                try:
                    out_polys.extend(poly.intersection(clip))
                except:
                    pass
        else:
            out_polys.append(poly)
    
#   # Treat GMT polygons
#   from ...misc.plot import map as Map
#   for poly in gmt_polys:
#       if isinstance(poly, Basemap) and poly.resolution is not None:
#           # We already have all the numeric polygons
#           gmt_m = poly
#       elif isinstance(poly, str):
#           # We must deduce polygons from gmt resolution
#           print 'getting gmt map for mask'
#           assert poly[0] in ['f', 'h', 'i', 'l', 'c'], \
#               "Resolution must be one of ['f', 'h', 'i', 'l', 'c'],  not "+poly[0]
#           # What kind of map?
#           if m is not None:
#               proj = m.projection
#           else:
#               proj = 'cyl'
#           # Ok, get it
#           gmt_m = Map(lon=lon, lat=lat, maponly=True, res=poly[0], projection=proj)
#       else:
#           continue
#       out_polys.extend([Polygon(N.array(p).transpose()) for p in gmt_m.coastpolygons])
        
    return out_polys
    

def _trans_(trans, xx, yy):
    #FIXME: see deg2map
    # Nothing to do
    if trans is False: return xx, yy
    # Check if we are sure to not have degrees
    xx = N.asarray(xx)
    yy = N.asarray(yy)
    xmin = xx.min()
    ymin = yy.min()
    xmax = xx.max()
    ymax = yy.max()
    if xmin >= -0.1 and ymin >= -0.1 and (xmax > 720. or ymax > 90.): 
        return xx, yy
    # We use a existing basemap instance
    if isinstance(trans, Basemap):
        return trans(xx, yy)
    # We create an instance
    lat_center = yy.mean()
    m = Basemap(lat_0=lat_center, lat_ts=lat_center, lon_0=xx.mean(), resolution=None,projection='merc',
        llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax)
    return m(xx, yy)
#   return d2m(xx, yy)

def d2m(x, y):
    return deg2m(x, y), deg2m(y)
    
def polygon_select(xx, yy, polys, zz=None, mask=False,):
    """Select unstructered points that are inside polygons
    
    - **xx**: Positions along X.
    - **yy**: Positions along Y.
    - **polys**: Polygons.
    - *zz*: Values at theses positions.
        
    :Outputs:
    
        - if ``mask is False``: selected ``x,y`` or ``x,y,z`` 
        - if ``mask is True``:
        
            - if ``zz is None``: the mask
            - else: all ``x,y,z`` where ``z`` is masked
    """

    # Get numeric axes
    xx = N.asarray(xx, 'd')
    yy = N.asarray(yy, 'd')

    # Get a proper polygon list
    polys = polygons(polys, lon=(xx.min(), xx.max()), lat=(yy.min(), yy.max()))
        
    # Create the mask
    pmask = N.zeros(xx.shape, '?')
    for ip, (x, y) in enumerate(zip(xx, yy)):
        for poly in polys:
            if Point((x,y)).within(poly):
                pmask[ip] = True
                break
    
    # Return mask or masked
    if mask:
        if zz is None:
            return pmask
        else:
            return xx, yy, N.ma.masked_where(pmask, zz, copy=0)
            
    # Return selection
    good = ~pmask
    del pmask
    xsel = xx[good]
    ysel = yy[good]
    ret = (xsel, ysel)
    if vv is not None:
        ret += (vv[good], )
    del good
    return ret


def convex_hull(xy, poly=False, method='delaunay'):
    """Get the envelop of cloud of points
    
    - **xy**: (x,y) or array of size (2,nxy)
    - *poly*: 
    
        - ``True``: Return as Polygon instance.
        - ``False``: Return two 1D arrays ``xpts,ypts``
        
        - *method*: 
        
            - ``"angles"``: Recursive scan of angles between points.
            - ``"delaunay"``: Use Delaunlay triangulation.
            
    """
    
    # Points
    xx, yy = M.get_xy(xy, m=False)
    xx = N.asarray(xx)
    yy = N.asarray(yy)
    if xx.ndim>1:
        xx = xx.ravel()
        yy = yy.ravel()
    
    # Various methods
    if method.startswith('delau'):
        
        # Delaunay
        from matplotlib.delaunay import Triangulation
        xy = uniq(N.asarray([xx, yy]).transpose())
        hull = Triangulation(xy[:, 0], xy[:, 1]).hull
        xe = xy[:, 0][hull]
        ye = xy[:, 1][hull]
        
    elif method == 'quarters':
        
        np = len(xx)
        
        # Most south point
        ip = N.argmin(yy)
        xe = [xx[ip]]
        ye = [yy[ip]]
        
        # SW
        while True:
            good = xx>xe[-1]
            if not good.any(): break
            ip = N.argmin(yy[good])
            xe.append(xx[ip])
            ye.append(yy[ip])
            
        # NW
        while True:
            good = yyx>ye[-1]
            if not good.any(): break
            ip = N.argmax(xx[good])
            xe.append(xx[ip])
            ye.append(yy[ip])
        
        pass
        #TODO: finish convex_hull with quaters
        
    elif method=='angles':
        
        # Angles
        np = len(xx)
        xx0 = N.resize(xx, (np, np))
        xx1 = N.resize(xx, (np, np)).transpose()
        dx = xx1-xx0 ; del xx0, xx1
        yy0 = N.resize(yy, (np, np))
        yy1 = N.resize(yy, (np, np)).transpose()
        dy = yy1-yy0 ; del yy0, yy1
        angles = N.arctan2(dx, dy) ; del dx, dy
        idx = N.arange(np)
        angles[idx, idx] = 10.
        
        # Most south point
        ip0 = ip = N.argmin(yy)
        xe = [xx[ip]]
        ye = [yy[ip]]
        
        # Recursive search
        ic = 0
        while True:
            ip = N.argmin(angles[ip])
            if ip == ip0: break
            xe.append(xx[ip])
            ye.append(yy[ip])
            ic += 1
            if ic > np: break
        xe = N.asarray(xe)
        ye = N.asarray(ye)
    
    else:
        raise NotImplementedError
    
    # Polygon or positions?
    if poly:
        return polygons([(xe, ye)], m=False)[0]
    return xe, ye


def uniq(data):
    """Remove duplicates in data
    
    :Example:
    
    >>> import numpy as N
    >>> a = N.arange(20.).reshape((10,2))
    >>> a[5] = a[1]
    >>> b = uniq(a)
    """
    if data.ndim>2:
        wdata = data.reshape((data.shape[0], data.size/data.shape[0]))
    else:
        wdata = data
    keep = N.ones(len(data), '?')
    for i, d in enumerate(data):
        if not keep[i]: continue
        bad = data == d
        if data.ndim>1:
            bad = bad.all(axis=-1)
        bad[i] = False
        keep = keep & ~bad
    return data[keep]
    

def envelop(*args, **kwargs):
    """Shortcut to :func:`convex_hull`"""
    return convex_hull(*args, **kwargs)

def grid_envelop(gg, centers=False, poly=True):
    """Return a polygon that encloses a grid
    
    - **gg**: A cdms grid or a tuple of (lat,lon)
    - *centers*: 
    
        - ``True``: The polygon is at the cell center. 
        - ``False``: The polygon is at the cell corners. 
        
    - *poly*: 
    
        - ``True``: Return as Polygon instance.
        - ``False``: Return two 1D arrays ``xpts,ypts``
    """
    # Axes
    x, y = M.get_xy(gg, num=True)
    
    # Rectangular grids
    gg = M.curv2rect(gg,mode=False)
    if M.isgrid(gg, curv=False):
        if centers:
            x0 = x.min() ;  x1 = x.max()
            y0 = y.min() ;  y1 = y.max()
        else:
            xb = M.bounds1d(x)
            yb = M.bounds1d(y)
            x0 = xb.min() ;  x1 = xb.max()
            y0 = yb.min() ;  y1 = yb.max()
        if poly:
            return Polygon(N.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], 'd'))
        return N.array([x0, x1]), N.array([y0, y1])
    
    # Get the 2D axes
    xx, yy = M.meshgrid(x, y)
    xxb, yyb = M.meshbounds(x, y)
    if centers:
        xxref, yyref = xx, yy
    else:
        xxref, yyref = xxb, yyb
    ny, nx = xx.shape
    
    # Get the limits
    xp = []
    yp = []
    for jslice, islice in (
        (0,  slice(None)), # bottom          
        (slice(1, -1),    -1),  # inner-left
        (-1, slice(None, None, -1)), # top
        (slice(-2, 0, -1), 0)): # inner right
            
        # Border
        xborder = xxref[jslice, islice].tolist()
        yborder = yyref[jslice, islice].tolist()        
            
        # Smart compression
        if isinstance(jslice, int): # Along X
            if jslice: # top
                jb0 = -2 ; jb1 = None
            else: # bottom
                jb0 = 0 ; jb1 = 2
            xpoly = xxb[jb0:jb1, islice][:, ::nx]
            ypoly = yyb[jb0:jb1, islice][:, ::nx]
        else: # Along Y
            if islice: # right
                ib0 = -2 ; ib1 = None
            else: # left
                ib0 = 0 ; ib1 = 2 
            xpoly = xxb[jslice, ib0:ib1][::ny-2]
            ypoly = yyb[jslice, ib0:ib1][::ny-2]
        xline = xx[jslice, islice] ; xline = xline[::len(xline)-1]
        yline = yy[jslice, islice] ; yline = yline[::len(yline)-1]
        poly1d = Polygon(N.asarray([[xpoly[0, 0], ypoly[0, 0]], 
            [xpoly[0, -1], ypoly[0, -1]], [xpoly[-1, -1], ypoly[-1, -1]], 
            [xpoly[-1, 0], ypoly[-1, 0]]], 'd'))
        line1d = LineString(N.asarray([[xline[0], yline[0]], [xline[1], yline[1]]], 'd'))
        if line1d.within(poly1d): # Compression
            xborder = xborder[::len(xborder)-1]
            yborder = yborder[::len(yborder)-1]
            
        # Add it
        xp += xborder
        yp += yborder
        
    # Return the polygon or arrays
    if not poly: return N.asarray(xp), N.asarray(yp)
    return Polygon(N.asarray(zip(xp, yp)))

def grid_envelop_mask(ggi, ggo, poly=True, **kwargs):
    """Create a mask on output grid from the bounds of an input grids:
    points of output grid that are not within these bounds of are masked.
    
    - **ggi**: Input grid or (lon,lat) or cdms variable.
    - **ggo**: Output grid or (lon,lat) or cdms variable.
    - Other keywords are passed to :func:`grid_envelop`
    
    :Returns:
        A mask on output grid.
    """
    
    if M.isgrid(ggi, curv=False): # Rectangular
        x, y = M.get_xy(ggi)
        xb = M.bounds1d(x)
        yb = M.bounds1d(y)
        ggo = M.get_grid(ggo)
        xxo, yyo = M.meshgrid(*M.get_xy(ggo))
        return (xxo<=xb.min())|(xxo>=xb.max())|(yyo<=yb.min())|(yyo>=yb.max())
    
    elif poly: # Using a polygon
    
        # Envelop of the input grid
        poly = grid_envelop(ggi, **kwargs)
        # Mask on output grid
        return ~polygon_mask(ggo, [poly], mode='inside')
        
    else: # Using regridding
        from regridding import regrid2d
        xxb, yyb = M.meshbounds(*M.meshbounds(*M.get_xy(ggi)))
        vari = MV2.asarray(xxb*0.)
        vari[[0, -1]] = 1.
        vari[:, [0, -1]] = 1.
        M.set_grid(vari, M.curv_grid(xxb, yyb))
        varo = regrid2d(vari, ggo, 'nearest', ext=True).filled(1.)
        return varo.astype('?')
        



def check_poly_islands(mask, polys, offsetmin=.85, offsetmax=1.5, dcell=2):
    """Check that islands as greater as a cell are taken into account
    
    - **mask**: The initial mask. A cdms variable with X (longitude) and Y (latitude), or a tuple (lon, lat, mask).
    - **polys**: Coastal polygons.
    - *offset*: Islands whose area is > 1-offset and < 1+offset % of the mean cell area are checked [default: .95]
    - *dcell*: dx and dy relative extension from the center of the cell in resolution units.
    
    """
    
    # Get mask and axes
    if isinstance(mask, tuple):
        if len(mask) == 3:
            xx, yy, mask = mask
        else:
            gg, mask = mask
            xx, yy = M.get_xy(gg)
    else:
        xx, yy = M.get_xy(mask)
    
    # Resolution
    dx = N.diff(xx).mean()
    dy = N.diff(yy).mean()
    cell_area = N.sqrt(dx*dy)
    
    # Loop on islands
    ni = 0
    for island_poly in polygons(polys):
        
        # Select islands
        island_area = island_poly.area()
        if island_area < cell_area*offsetmin or island_area > cell_area*offsetmax: continue
        
        # Define cells around
        xisland = island_poly.get_data()[:, 0]
        yisland = island_poly.get_data()[:, 1]
        xcisland = (xisland.min()+xisland.max())/2.
        ycisland = (yisland.min()+yisland.max())/2.
        imin, imax, k = xx.mapInterval((xcisland-dx*dcell*1.1, xcisland+dx*dcell*1.1))
        jmin, jmax, k = xx.mapInterval((ycisland-dy*dcell*1.1, ycisland+dy*dcell*1.1))
        amask = mask[jmin:jmax, imin:imax]
        xxb, yyb = _trans_(trans, bounds2d(amask.getAxis(-1)),  bounds2d(amask.getAxis(-2)))
        
        # Get intersection areas
        areas = N.zeros(amask.shape, 'f')
        for j in xrange(amask.shape[-2]):
            for i in xrange(amask.shape[-1]):
                cell_poly = Polygon(N.array([xxb[j, i, :].ravel(), yyb[j, i, :].ravel()]).transpose())
                if cell_poly.intersects(island_poly):
                    areas[j, i] = cell_poly.intersection(island_poly).area()
                    
        # It should be masked where area is max
        if not amask.flat[N.argmax(areas)]: ni += 1
        amask.flat[N.argmax(areas)] = True
        mask[jmin:jmax, imin:imax] = amask
        
    print 'Masked %i islands' % ni
    return mask
    
def check_poly_straits(mask, polys, dcell=2, threshold=.75):
    """Check that straits are opened by counting the number of polygon intersection in each cell and the area of water
    
    - **mask**: The initial mask. A cdms variable with X (longitude) and Y (latitude).
    - **polys**: Coastal polygons.
    - *dcell*: dx and dy relative extension from the center of the cell in resolution units.
    - *threshold*: Maximal fraction of cell area that must be covered of land (> .5) [default: .25]
    """
    
    # Get mask and axes
    if isinstance(mask, tuple):
        if len(mask) == 3:
            xx, yy, mask = mask
        else:
            gg, mask = mask
            xx, yy = M.get_xy(gg)
    else:
        xx, yy = M.get_xy(mask)
        
    # Get cell bounds
    xxb, yyb = bounds2d(xx),  bounds2d(yy)
    
    # Loop on coastal points
    ns = 0
    for i, j in get_coastal_indices(mask):
        
        # Cell polygon
        cell_poly = Polygon(N.array(xxb[j, i, :], yyb[j, i, :]).transpose())
        cell_area = cell_poly.area()
        
        # Count polygon intersections
        npoly = 0
        area = 0.
        for poly in polys:
            if poly.intersects(cell_poly):
                intersections = poly.intersection(cell_poly)
                npoly += len(intersections)
                for intersection in intersections:
                    area += intersection.area()
                    if area > threshold: break
                else:
                    continue
            break
        else:
            # Several polygons = strait, so open
            mask[j, i] = npoly < 2
            
    print 'Opened %i straits' %ns
    return mask
    

def t2uvmasks(tmask, getu=True, getv=True):
    """Compoute masks at U and V points from a mask at T points on a C grid (True is land)
    
    - **tmask**: A 2D mask.
    - *getu*: Get the mask at U-points
    - *getv*: Get the mask at V-points
    
    Return umask,vmask OR umask OR vmask depending on getu and getv.
    """
    if not getu and not getv: return
    if getu:
        umask = N.array(tmask) # copy
        umask[..., :, :-1] = tmask[..., :, :-1] | tmask[..., :, 1:]
        if not getv: return umask
    vmask = N.array(tmask)
    vmask[..., :-1, :] = tmask[..., :-1, :] | tmask[..., 1:, :]
    if not getu: return vmask
    return umask, vmask

def mask2d(mask, method='and'):
    """Convert a 3(or more)-D mask to 2D by performing compression on first axes
    
    .. note::
    
        Mask is False on ocean
    
    - **mask**: At least 2D mask
    - *method*: Compression method
    
        - ``"and"``: Only one point must be unmask to get the compressed dimension unmasked.
        - ``"or"``: All points must be unmask to get the compressed dimension unmasked.
    """
    assert mask.ndim > 1, 'Mask must be at leat 2D'
    if mask.ndim == 2: return mask
    dtype = mask.dtype
    if method.lower() == 'or':
        func = N.logical_or
    else:
        func = N.logical_and
    while mask.ndim != 2:
        mask = func.reduce(mask)
    return mask.astype(dtype)
    
    
def _old_rsamp_(x, y, r, z=None, rmean=.7, proj=False, getmask=False):
    """Undersample data points using a radius of proximity
    
    - **x**: 1D x array.
    - **y**: 1D Y array.
    - **r**: Radius of proximity.
    - *z*: 1D Z array.
    - *proj*: Geographic projection instance to convert coordinates.
    - *method*:
    
        - ``mean``: Average neighbouring points.
        - ``pickup``: Pickup only the current point within the circle of proximity.
    """
    # Need numpy arrays
    x = N.asarray(x)
    y = N.asarray(y)
    if getmask:
        z = None
    if z is not None:
        z = N.array(z)
    n = x.shape[0]
    dst = x.copy()
    dst[:] = 0.
    close = N.zeros(n, '?')
    good = N.ones(n, '?')
    rmean = N.clip(rmean, 0., 1.)
    
    # Projection
    if proj is True:
        proj = get_proj((x,y))
    if callable(proj):
        x, y = proj(x, y)
        
    # Loop on valid points
    for i0 in xrange(n):
        if not good[i0]: continue
        # Search for neighbours
        dst[:] = N.sqrt((x-x[i0])**2+(y-y[i0])**2)
        close[:] = dst < r
        # Average over neighbours
        if z is not None and rmean>0.:
            zclose = z[close]
            z[i0] = zclose[dst[close] < r*rmean].mean()
            del zclose
        # Neighbours no longer valid
        good &= ~close
        good[i0] = True
    del dst, close
    # Only the mask
    if getmask:
        return good
    # Select valid points
    ret = (x[good], y[good])
    if z is not None:
        ret += z[good], 
    del good
    return ret

def rsamp(x, y, r, z=None, rmean=.7, proj=False, rblock=3, getmask=False):
    """Undersample data points using a radius of proximity
    
    - **x**: 1D x array.
    - **y**: 1D Y array.
    - **r**: Radius of proximity.
    - *z*: 1D Z array.
    - *proj*: Geographic projection instance to convert coordinates.
    - *rmean*: Radius of averaging relative to ``r`` (``0<rmean<1``)
    - *rblock*: Size of blocks relative to ``r`` for computation by blocks.
    """
    # Need numpy arrays
    x = N.asarray(x)
    y = N.asarray(y)
    if getmask:
        z = None
    if z is not None:
        z = N.array(z)
    n = x.shape[0]
    good = N.ones(n, '?')
    rmean = N.clip(rmean, 0., 1.)
    rblock = int(max(1, rblock))
    
    # Projection
    if proj is True:
        proj = get_proj((x,y))
    if callable(proj):
        x, y = proj(x, y)

    # Setup blocks
    xmin = x.min()
    ymin = y.min()
    xmax = x.max()
    ymax = y.max()
    dx = xmax-xmin ; dy = ymax-ymin 
    xmin -= .001*dx ; xmax += .001*dx
    ymin -= .001*dy ; ymax += .001*dy
    nxb = max(1, int((xmax-xmin)/(r*(rblock-1))))
    nyb = max(1, int((ymax-ymin)/(r*(rblock-1))))
    rblock *= r
    xgood = N.zeros(n, '?')
    block = N.zeros(n, '?')
    
    # Loop on x bands
    x0 = xmin+r-rblock
    for ixb in xrange(nxb):
        x0 += rblock-r
        if ixb == nxb-1:
            x1 = xmax
        else:
            x1 = x0 + rblock
        xgood[:] = (x>=x0)&(x<=x1)
        if not xgood.any(): continue
        
        # Loop on y bands
        y0 = ymin+r-rblock
        for iyb in xrange(nyb):
            y0 += rblock-r
            if iyb == nyb-1:
                y1 = ymax
            else:
                y1 = y0 + rblock
                
            # Select the block
            block[:] = xgood&(y>=y0)&(y<=y1)
            xx = x[block]
            nn = len(xx)
            if not nn: continue
            yy = y[block]
            if z is not None:
                zz = z[block]
            gg = good[block]
            
            # Loop on valid points
            for i0 in xrange(nn):
                if not gg[i0]: continue
                # Search for neighbours
                dst = N.sqrt((xx-xx[i0])**2+(yy-yy[i0])**2)
                close = dst < r
                # Average over neighbours
                if z is not None and rmean>0.:
                    zclose = zz[close]
                    zz[i0] = zclose[dst[close] < r*rmean].mean()
                    del zclose
                # Neighbours no longer valid
                gg &= ~close
                gg[i0] = True
                # del
                del dst, close
                
            # Store results
            good[block] = gg
            if z is not None:
                z[block] = zz
                del zz
            del xx, yy, gg
    # Only the mask
    if getmask:
        return good
    # Select valid points
    ret = (x[good], y[good])
    if proj:
        ret = proj(ret[0], ret[1], inverse=True)
    if z is not None:
        ret += z[good], 
    del good
    return ret


def zcompress(z, *xy, **kwargs):
    """Compress 1D arrays according to the mask of the first one
    
    - **z**: Reference (masked) array
    - *xy*: Additional arrays to be compressed
    - *numpify*: Force convertion to numpy array
    
    :Example:
    
    >>> z, x, y = zcompress(z, x, y, numpify=True)
    """
    if hasattr(z, 'mask') and z.mask is not N.ma.nomask:
        good = ~z.mask
        ret = ()
        npf = kwargs.get('numpify', False)
        for var in (z, )+xy:
            # Compress
            if not cdms2.isVariable(var):
                newvar = var[good]
            else:
                newvar = cdms2.createVariable(var.asma()[good])
                cp_atts(var, newvar, id=True)
            # Numpify
            if npf and hasattr(var, 'filled'):
                ret += newvar.filled(), 
                del newvar
            else:
                ret += newvar, 
        return ret
    ret = z, 
    ret += xy
    return ret
   
def resol_mask(grid, res=None, xres=None, yres=None, xrelres=None, yrelres=None, 
    relres=None, scaler=None, compact=False, relmargin=.05):
    """Create a mask based on resolution criteria for undersampling 2D data
    
    :Params:
        - **grid**: A cdms array with a grid, a cdms grid or a tuple of axes.
        - **res**, optional: Horizontal resolution of arrows 
            (in both directions) for undersampling [default: ``None``]
            
            .. warning::
    
                Use of ``res`` makes the supposition that X and y units are consistent.
        
        - **xres**, optional: Same along X [default: ``res``]
        - **yres**, optional: Same along Y [default: ``res``]
        - **relres**, optional: Relative resolution 
            (in both directions). 
            
            - If > 0, = ``median(res)*relres``. 
            - If < -1, =``min(res)*abs(relres)``. 
            - If < 0 and > -1, =max(res)*abs(relres)
            
        - **xrelres**, optional: Same along X [default: ``relres``]
        - **yrelres**, optional: Same along Y [default: ``relres``]
        - **scaler**, optional: A callable object that transform X and Y coordinates.
          Transformed coordinates are used instead of original coordinates
          if resolutions are negatives.
          A typical scaler is for example a geographic projection that convert degrees
          to meters (like a :class:`~mpl_toolkits.basemap.Basemap` instance).
          
        - **compact**, optional: If no unsersampling is efective, returns ``False``.
          
    .. note:: A resolution value set to ``False`` implies no undersampling.
          
    :Example:
    
        >>> mask = resol_mask((x2d, y2d), relres=2.5) # sampling=2.5
        >>> mask = resol_mask(var.getGrid(), xres=2., scaler=mymap, yres=-100.) # 100m
          
    .. note::
    
        It is usually better to regrid with a convenient method 
        instead of unsersample because of aliasing.
        
    """
    # Get numeric axes
    xx, yy = M.get_xy(grid, num=True, mesh=True, m=False)
    
    # Scaler
    if callable(scaler):
        xxs, yys = scaler(xx, yy)
    else:
        xxs, yys = xx, yy
    
    # Guess reference resolution
    if res is not None:
        if xres is None: xres = res
        if yres is None: yres = res
    if relres is not None:
        if xrelres is None: xrelres = relres
        if yrelres is None: yrelres = relres        
    if xres is False: xres = None
    if yres is False: yres = None
    if xrelres is False: xrelres = None
    if yrelres is False: yrelres = None
    if (xrelres is not None and xres is None) or (yrelres is not None and yres is None):
        xmode = 'median'
        ymode = 'median'
        if xrelres is not None:
            if xrelres < -1:
                xmode = 'min'
            elif xrelres < 0:
                xmode = 'max'
        if yrelres is not None:
            if yrelres < -1:
                ymode = 'min'
            elif yrelres < 0:
                ymode = 'max'
        xr = M.resol(xx, axis=1, mode=xmode) if xres is None else xres
        yr = M.resol(yy, axis=0, mode=ymode) if yres is None else yres
        if N.abs(N.log10(xres/yres))<2: # similar coordinates
            xr, yr = M.resol(xxs, yys, mode=(xmode, ymode))
            if xres is None: xres = -xr
            if yres is None: yres = -yr
        else:
            if xres is None: xres = xr
            if yres is None: yres = yr
        
    # Resolution along X
    if xres is not None: # From absolute resolution
        if xrelres is None: xrelres = 1.
        x2d = xx if xres > 0 else xxs
        xmask = N.zeros(xx.shape, '?')
        baserelres = M.resol(x2d, axis=1, mode='loc')/N.abs(xres*N.abs(xrelres))
        newrelres = baserelres.copy()
        level = 1
        while level < 100:
            bad = newrelres < 1-relmargin
            if not bad.any(): break
            levelmask = bad.copy()
            levelmask[:, level/2::level+1] = False
            xmask = N.where(bad, levelmask, xmask)
            newrelres[bad] = baserelres[bad]*(level+1)
            del bad, levelmask
            level += 1
        del baserelres, newrelres
    else:
        xmask = False
            
    # Resolution along Y
    if yres is not None: # From absolute resolution
        if yrelres is None: yrelres = 1.
        y2d = yy if yres > 0 else yys
        ymask = N.zeros(yy.shape, '?')
        baserelres = M.resol(y2d, axis=0, mode='loc')/N.abs(yres*N.abs(yrelres))
        newrelres = baserelres.copy()
        level = 1
        while level < 100:
            bad = newrelres < 1-relmargin
            if not bad.any(): break
            levelmask = bad.copy()
            levelmask[level/2::level+1, :] = False
            ymask = N.where(bad, levelmask, ymask)
            newrelres[bad] = baserelres[bad]*(level+1)
            del bad, levelmask
            level += 1
        del baserelres, newrelres
    else:
        ymask = False
    
    return xmask | ymask
 



######################################################################
######################################################################
import misc as M
from ...misc.axes import islon, islat
from ...misc.phys.units import deg2m
from ...misc.io import Shapes
from basemap import GSHHS_BM, get_proj
from ...misc.misc import cp_atts
