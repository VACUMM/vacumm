"""Test the fortran function :f:func:`curv2rel`"""
from vcmq import N, P, code_file_name, P, os, meshbounds
from vacumm.misc.grid._interp_ import curv2rel


# Input grid
x0, y0 = 0., 2.
nxi = 3
nyi = 3
dxi = (2., 2.)
dyi = (-2., 2.)
xxi = N.zeros((nyi, nxi))
xxi[0] = x0 + N.arange(nxi) * dxi[0]
for j in range(1, nyi):
    xxi[j] = xxi[j-1] + dyi[0]
yyi = N.zeros((nyi, nxi))
yyi[:, 0] = y0 + N.arange(nyi) * dyi[1]
for j in range(1, nyi):
    yyi[:, j] = yyi[:, j-1] + dxi[1]
xxbi, yybi = meshbounds(xxi, yyi)
relpos2index = lambda fi, fj: fj * nyi + fi
ii, jj = N.meshgrid(N.arange(nxi)+.5, N.arange(nyi)+.5)
iib, jjb = meshbounds(ii, jj)
zzi = jj * nyi + ii #relpos2index(ii, jj)
zzbi = jjb * nyi + iib # ny+1? relpos2index(iib, jjb)

# Input random points
N.random.seed(0)
np = 100
xxo = N.random.random(np)*(xxbi.max()-xxbi.min()) + xxbi.min()
yyo = N.random.random(np)*(yybi.max()-yybi.min()) + yybi.min()
xxo = N.concatenate((xxo, xxi.ravel()))
yyo = N.concatenate((yyo, yyi.ravel()))
np = xxo.size

# Convert to relative indices
pp, qq = curv2rel(xxbi, yybi, xxo, yyo)
bad = pp<0
zzo = (qq-1) * nyi + pp-1 #relpos2index(pp-1, qq-1)
zzo[bad] = -1

# Plot
vmin = zzbi.min()
vmax = zzbi.max()
P.pcolor(xxbi, yybi, zzi, vmin=vmin, vmax=vmax)
P.scatter(xxo, yyo, c=zzo, vmin=vmin, vmax=vmax, s=80)
P.grid()
P.title('curv2rel')
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
P.savefig(figfile)
P.close()
