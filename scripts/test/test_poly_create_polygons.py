"""Test the :func:`~vacumm.misc.poly.create_polygons` function"""
from vcmq import N, create_polygons, create_polygon, plot_polygon, create_grid
from _geoslib import Polygon

# Data
xx = N.array([0., 5., 4., 1.])
yy = N.array([0., 0., 2., 2.])
clip = [-3, -3, 0, 0]
xaxis = yaxis = N.linspace(-2., 2., 20.)

def plot_polygons(polys, **kwargs):
    for p in polys:
        plot_polygon(p, **kwargs)

# Single known argument
pp0 = create_polygon(Polygon(N.array([xx, yy]).T))

# Classic data
pp0 = create_polygons([N.array([xx, yy]), N.array([xx, yy]).T+6.])

# Classic with projection
proj = lambda x, y: (x*1.5, y*1.5)
pp1 = create_polygons([N.array([xx, yy])])

# Classic with clipping
pp2 = create_polygons([-N.array([xx, yy])], clip=clip)

## From grid
#pp3 = create_polygons([create_grid(xaxis, yaxis)]) # no longer integrated

