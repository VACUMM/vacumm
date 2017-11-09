"""Test the :func:`~vacumm.misc.grid.io.Shapes` class"""
from vcmq import N, P, data_sample
from vacumm.misc.io import Shapes
result = []

# Shapefile
shpfile = data_sample("ne_110m_land/ne_110m_land")

# Basic
S = Shapes(shpfile)
result.append(('assertEqual', [len(S), 127]))
S.plot(title='Basic', show=False, m_fig='new')

# Min area
S = Shapes(shpfile, min_area=1000)
result.append(('assertEqual', [len(S), 3]))
S.plot(title='Min area', show=False, m_fig='new')

# Projection
S = Shapes(shpfile, proj='merc')
result.append(('assertGreater', [S[1].area(), 1e12]))
S.plot(title='Projected', show=False, m_fig='new')

# Clips
S = Shapes(shpfile, clip=[-10, 42, 10, 51.])
result.append(('assertEqual', [len(S), 3]))
S.plot(title='Clipped', show=False, m_fig='new')
S = Shapes(shpfile, clip=[-10, 42, 10, 51.], proj=True)
result.append(('assertEqual', [len(S), 3]))
S.plot(title='Clipped+projected', show=False, m_fig='new')

P.show()

