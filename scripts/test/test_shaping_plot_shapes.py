"""Test the :func:`vacumm.misc.shaping.plot_shapes` function"""
from vcmqm import (N, P, create_shape, create_shapes, plot_shapes)

# %% Create shapes
# Points
xmin, xmax, ymin, ymax = -10, 43, 10, 65
xp = N.random.uniform(xmin, xmax, 10)
yp = N.random.uniform(ymin, ymax, 10)
points = create_shapes(N.array([xp, yp]).T, 'point')
# Lines
lines = [create_shape([[xmin, ymin], [xmax, ymin],
                       [xmin, ymax], [xmin, ymin]], 'linestring')]
# Polygons
polys = [create_shape([[xmin, ymin], [xmax, ymin], [xmax, ymax]], 'polygon')]

# %% Plots
P.figure(figsize=(6, 8))
P.subplot(311)
plot_shapes(points, color='tab:blue')
P.subplot(322)
plot_shapes(lines, color='tab:green', linewidth=3)
P.subplot(333)
plot_shapes(polys, color='.5', linewidth=0, points=True,
            points_c='tab:red', s=100, points_zorder=10)
