"""Test the :func:`~vacumm.misc.poly.create_polygons` function"""
from vcmq import (N, create_polygons, create_polygon,
                  clip_shapes)
from _geoslib import Polygon

# %% Data
xx = N.array([0., 5., 4., 1.])
yy = N.array([0., 0., 2., 2.])
clip = [-3, -3, 0, 0]
xaxis = yaxis = N.linspace(-2., 2., 20.)

# %% Single known argument
pp0 = create_polygon(Polygon(N.array([xx, yy]).T))

# %% Classic data
pp0 = create_polygons([N.array([xx, yy]), N.array([xx, yy]).T+6.])

# %% Classic with projection
proj = lambda x, y: (x*1.5, y*1.5)
pp1 = create_polygons([N.array([xx, yy])])

# %% Classic with clipping
pp2 = create_polygons([-N.array([xx, yy])])
pp2 = clip_shapes(pp2, clip)
