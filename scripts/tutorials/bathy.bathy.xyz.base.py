# %% Imports
from __future__ import print_function
import numpy as N, pylab as P
from vcmq import XYZBathy, create_grid, resol, map2

# %% Fake random bathy
n = 1000
xc = -8.
yc = 48.
xr = 5.
yr = 2.
noise = N.random.random(n)
a = N.random.random(n)*N.pi*2
r = N.random.random(n)
x = xc + r*xr*N.cos(a)
y = yc + r*yr*N.sin(a)
bathy = N.asarray([x, y, N.exp(-(x-xc)**2/xr**2 -
                               (y-yc)**2/yr**2)*30.+noise]).transpose()
fbathy = __file__[:-2]+'xyz'
N.savetxt(fbathy, bathy)

# %% Init a XYZBathy with some undersampling
xyz = XYZBathy(fbathy,  long_name='My XYZ', rsamp=0.2)

# %% Add a triangular selection area
xyz.select([[-12., 46], [0, 48], [-8, 50]])

# %% Add an exclusion area ([xmin,ymin,xmax,ymax])
xyz.exclude([-9, 47., -7, 48.25])

# %% Infos
print(xyz)

# %% Get values
x = xyz.x
y = xyz.y
z = xyz.z

# %% Extents
# - taking into account exclusions
print('Limits:', xyz.xmin,  xyz.xmax, xyz.ymin,  xyz.ymax)
# - raw data
print('Rax X min:', xyz.get_xmin(mask=False),   xyz.get_x(mask=False).min())

# %% Mean resolution
print('Resolution in degrees:', xyz.resol(deg=True))
print('Resolution in meters:', xyz.resol())

# %% Automatic definition of a rectangular grid
grid_auto = xyz.grid
print('Resolution of automatic grid:', resol(grid_auto))

# Interpolation on automatic grid
print('Interpolation auto')
gridded_auto = xyz.togrid(mask='l')
#  equivalent to:
#  >>> gridded_auto = xyz.togrid(xyz.grid, mask='l')

# %% Interpolation sur grille manuelle
print('Interpolation on user grid, masking and extraction')
# - define the gird
grid_manual = create_grid((-10, -3, .2), (47, 49, .2))
# - interpolation
gridded_manual = xyz.togrid(grid_manual, mask='l')
# - extraction
xyz_up = xyz.clip(zone=(None, None, None, 48.3), margin=2)
# - margin : relative margin in resolution units
#   -> ici : ymax = 48.3 + xyz.resol()[1]*2

# %% Dump to file
print('Save')
prefix = __file__[:-2]+'up'
xyz_up.save(prefix+'.xyz')  # ascii
xyz_up.save(prefix+'.nc')   # netcdf/grd

# %% Plots
print('Plots')
# - init
P.figure(figsize=(4.5, 8))
P.rc('font', size=8)
P.subplots_adjust(top=.95, hspace=.25, left=.1, bottom=.05, right=.98)
m = map2(lon=(xc-xr, xc+xr), lat=(yc-yr, yc+yr), proj='merc',
         subplot=311, autoresize=0, resolution='l', show=False,
         ticklabel_size=9, xhide=False)
kwplot = dict(vmin=xyz.get_zmin(False), vmax=xyz.get_zmax(False),
              m=m, show=False, colorbar=False)
xyz.plot(size=10, mode='both', masked_alpha=.1, **kwplot)
kwplot.update(autoresize=0, ticklabel_size=9)
map2(gridded_manual, subplot=312, xhide=False, title='Sur grille manuelle',
     **kwplot)
map2(gridded_auto, subplot=313, title='Sur grille auto',
     close=False, **kwplot)
