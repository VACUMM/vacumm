"""Test :func:`~vacumm.misc.grid.misc.isrect`"""

# %% Imports
from vcmq import isrect, create_grid, rotate_grid

# %% Trully rectangular
rgrid = create_grid((0., 25, .25), (0., 10, .25), 'rect')
assert isrect(rgrid)

# %% Rectangular 2D
crgrid = create_grid((0., 25, .25), (0., 10, .25), 'curv')
assert isrect(crgrid)

# %% Curvilinear
cgrid = rotate_grid(crgrid, 30.)
assert not isrect(cgrid)
