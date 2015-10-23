"""Test the :func:`~vacumm.misc.grid.io.Shapes` class"""
from vcmq import N, P, code_file_name, data_sample
from vacumm.misc.io import Shapes
result = []

# File names
shpfile = data_sample("ne_110m_land/ne_110m_land")
figfile = code_file_name(ext=False)+'_%i.png'

# Basic
S = Shapes(shpfile)
result.append(('assertEqual', [len(S), 127]))
S.plot(title='Basic', show=False)
P.savefig(figfile%0);P.close()

# Min area
S = Shapes(shpfile, min_area=1000)
result.append(('assertEqual', [len(S), 3]))
S.plot(title='Min area', show=False)
P.savefig(figfile%1);P.close()

# Projection
S = Shapes(shpfile, proj='merc')
result.append(('assertGreater', [S[1].area(), 1e12]))
S.plot(title='Projected', show=False)
P.savefig(figfile%2);P.close()

# Clips
S = Shapes(shpfile, clip=[-10, 42, 10, 51.])
result.append(('assertEqual', [len(S), 3]))
S.plot(title='Clipped', show=False)
P.savefig(figfile%3);P.close()
S = Shapes(shpfile, clip=[-10, 42, 10, 51.], proj=True)
result.append(('assertEqual', [len(S), 3]))
S.plot(title='Clipped+projected', show=False)
P.savefig(figfile%4);P.close()

