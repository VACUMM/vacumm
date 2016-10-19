"""Test the fortran function :f:func:`nearest2d`"""
from vcmq import N, P, code_file_name, P, os, rotate_grid, add_grid, meshbounds
from vacumm.misc.grid._interp_ import nearest2d


# Input grid
gridi = rotate_grid((N.arange(5), N.arange(4)), 30)
xxi = gridi.getLongitude()[:].filled()
yyi = gridi.getLatitude()[:].filled()
vari = N.resize(yyi, (20, )+ yyi.shape)
nb = 10
xxbi, yybi = meshbounds(xxi, yyi)

# Output grid
grido = rotate_grid((N.linspace(0, 6, 50)-1, N.linspace(0, 4, 35)+1.), -20)
xxo = grido.getLongitude()[:].filled()
yyo = grido.getLatitude()[:].filled()
xxbo, yybo = meshbounds(xxo, yyo)

# Nearest
varo = nearest2d(vari, xxi, yyi, xxo, yyo, nb)

# Plot
vmin = varo.min()
vmax = varo.max()
P.figure(figsize=(8, 4))
P.subplot(121, aspect=1)
P.pcolor(xxbi, yybi, vari[0], vmin=vmin, vmax=vmax)
add_grid(grido)
P.title('original')
P.subplot(122, aspect=1)
P.pcolor(xxbo, yybo, varo[0], vmin=vmin, vmax=vmax)
add_grid(gridi)
P.title('nearest2d')
P.axis('image')
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
P.savefig(figfile)
P.show()
P.close()
