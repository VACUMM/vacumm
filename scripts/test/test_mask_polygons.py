"""Test the :func:`~vacumm.misc.grid.masking.polygons` function"""
from vcmq import N, P, polygons, create_polygon, plot_polygon, create_grid
from _geoslib import Polygon

# Data
xx = N.array([0., 5., 4., 1.])
yy = N.array([0., 0., 2., 2.])
clip = [-3, -3, 0, 0]
xaxis = yaxis = N.linspace(-2., 2., 20.)

def plot_polygons(polys, **kwargs):
    for p in polys:
        plot_polygon(p, **kwargs)

## Single known argument
#pp0 = polygons(Polygon(N.array([xx, yy]).T))
#
## Classic data
#pp0 = polygons([N.array([xx, yy]), N.array([xx, yy]).T+6.])
#
## Classic with projection
#proj = lambda x, y: (x*1.5, y*1.5)
#pp1 = polygons([N.array([xx, yy])])
#
## Classic with clipping
#pp2 = polygons([-N.array([xx, yy])], clip=clip)

# From grid
pp3 = polygons([create_grid(xaxis, yaxis)])

