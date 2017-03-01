"""Test the fortran function :f:func:`curv2rel_single`"""
from vcmq import N, P, code_file_name, P, os, meshbounds
from vacumm.misc.grid._interp_ import curv2rel_single


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
relpos2index = lambda fi, fj, nyi: fj * nyi + fi
ii, jj = N.meshgrid(N.arange(nxi)+.5, N.arange(nyi)+.5)
iib, jjb = meshbounds(ii, jj)
zzi = relpos2index(ii, jj, nyi)
zzbi = relpos2index(iib, jjb, nyi)

# Input random points
N.random.seed(0)
np = 100
xxo = N.random.random(np)*(xxbi.max()-xxbi.min()) + xxbi.min()
yyo = N.random.random(np)*(yybi.max()-yybi.min()) + yybi.min()
xxo = N.concatenate((xxo, xxi.ravel()))
yyo = N.concatenate((yyo, yyi.ravel()))
np = xxo.size

# Convert to relative indices
zzo = N.zeros(np)-1
for i, (xo, yo) in enumerate(zip(xxo, yyo)):
    p, q = curv2rel_single(xxbi, yybi, xo, yo)
#    print xo, yo, '|', p, q, (q-1) * nxi + p - 1
    if q>=0:
        zzo[i] = relpos2index(p-1, q-1, nyi)

# Plot
vmin = zzbi.min()
vmax = zzbi.max()
P.pcolor(xxbi, yybi, zzi, vmin=vmin, vmax=vmax)
P.scatter(xxo, yyo, c=zzo, vmin=vmin, vmax=vmax, s=80)
P.grid()
P.title('curv2rel_single')
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
P.savefig(figfile)
P.close()
