"""Test the :func:`~vacumm.misc.grid.masking.create_polygon` function"""
from vcmq import N, P, create_polygon, plot_polygon, code_file_name
from _geoslib import Polygon

# Data
xx = N.array([0., 5., 4., 1.])
yy = N.array([0., 0., 2., 2.])

# From Polygon
P.figure(figsize=(5, 5))
p0 = create_polygon(Polygon(N.array([xx, yy]).T))
plot_polygon(p0, color='b', label='From poly')

# From data
p1 = create_polygon(N.array([xx, yy]).T-.25)
plot_polygon(p1, color='r', label='From data')

# From transposed data 
p2 = create_polygon(N.array([xx, yy])-0.5)
plot_polygon(p2, color='g', label='From transposed data')

# From min/max
p3 = create_polygon([xx.min(), yy.min(), xx.max(), yy.max()])
plot_polygon(p3, color='m', label='From min/max')

# With projection
proj = lambda x, y: (x*1.5, y*1.5)
p4 = create_polygon(p0, proj=proj)
plot_polygon(p4, color='cyan', label='With projection')

# Save
P.title('create_polygon')
P.legend(loc='best')
P.savefig(code_file_name(ext='png'))
P.show()
P.close()
