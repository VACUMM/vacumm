# Init
from vacumm.misc.io import write_snx
from _geoslib import Point, LineString, Polygon
import numpy as N
fbase = __file__[:-2]

# Points
fpoints = fbase+'points.snx'
# - tuples de points
write_snx([(0., 0., 0.), (1., 1., 1.)], fpoints)
# - marche aussi par groupe
write_snx([[(2., )*3, (3, )*3]], fpoints, mode='a') # append
# - depuis des tableaux
x = N.arange(4.)
y = N.arange(4.)
z = N.arange(4.)
bloc_xyz = N.array([x, y, z]).transpose() # bloc_xyz[i] = point i
write_snx(bloc_xyz, fpoints, mode='a')
# - deux blocs
bloc_xyz2 = bloc_xyz+10.
write_snx([bloc_xyz, bloc_xyz2], fpoints, type='points', mode='a')
#   note : il faut preciser le type dans ce cas
# - depuis des "Point" de _geoslib
write_snx([Point((0., 0.)), Point((100., 100.))], fpoints, mode='a', z = 99)

# Lignes
flines= fbase+'lines.snx'
# - blocs
write_snx(bloc_xyz, flines, type='lines') # on doit preciser le type
write_snx([bloc_xyz, bloc_xyz2], flines, mode='a')
# - "LineString" de _geoslib
write_snx(LineString(bloc_xyz[:, :2]+100.), flines, mode='a', z=999)
write_snx([LineString(bloc_xyz[:, :2]+100.), LineString(bloc_xyz2[:, 2]+100.)], flines, mode='a')

# Polygones : comme les lignes, mais fermes
fpolys= fbase+'polys.snx'
# - blocs
write_snx(bloc_xyz, fpolys, type='polygons') # on doit preciser le type
write_snx([bloc_xyz, bloc_xyz2], fpolys, mode='a', type='polygon')
# - on referme presque la ligne => detection AUTO comme polygone
bloc_xyz3 = N.concatenate((bloc_xyz, bloc_xyz[:1]+1.))
write_snx([bloc_xyz3], fpolys, mode='a')
# - "Polygon" de _geoslib
write_snx(Polygon(bloc_xyz3[:, :2]+100.), flines, mode='a')

# Via un descripteur de fichier
f = open(fbase+'mix.snx', 'w')
write_snx([bloc_xyz, bloc_xyz2], f, type='point', close=False)
write_snx([bloc_xyz, bloc_xyz2], f, type='line', close=False)
write_snx([bloc_xyz, bloc_xyz2], f, type='polygon', close=False)
f.close()

# Ecriture auto dans des fichiers separes
write_snx([bloc_xyz, bloc_xyz2], fbase+'split%i.snx', type='polygon', close=False)


print 'Done'
