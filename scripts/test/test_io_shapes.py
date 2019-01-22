"""Test the :class:`~vacumm.misc.sdata.Shapes` class"""
from vcmqm import P, data_sample, Shapes
import pyproj

# %% Shapefile
shpfile = data_sample("ne_110m_land/ne_110m_land")

# %% Basic
S = Shapes(shpfile)
assert len(S) == 127
S.plot(title='Basic', show=False, fig=P.figure())

# %% Projection
x, y = S.get_xy(split=True)
proj = pyproj.Proj({'proj':'moll', 'lon_0':0.})
S = S.transform(proj)
assert S[0].area() > 1e12
S.plot(title='Projected', show=False, fig=P.figure())

# %% Clips
S = Shapes(shpfile, clip=[-10, 42, 10, 51.])
assert len(S) == 3
S.plot(title='Clipped', show=False, fig=P.figure())
S = Shapes(shpfile, clip=[-10, 42, 10, 51.]).transform(proj)
assert len(S) == 3
S.plot(title='Clipped+projected', show=False, fig=P.figure())
