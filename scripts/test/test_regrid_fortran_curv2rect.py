"""Test the fortran function :f:func:`curv2rect`"""
from vcmq import N, P, code_file_name, P, os
from vacumm.misc.grid._interp_ import curv2rect

# Input
x1, y1 = 0., 0.
x2, y2 = 3., 1.
x3, y3 = 2., 4.
x4, y4 = -1., 2.

# Format and convert
xx, yy = N.meshgrid(N.arange(-2, 4, 0.25), N.arange(-1, 5, 0.25))
nxy = xx.shape
xx.shape = -1
yy.shape = -1
pp, qq = [], []
for x, y in zip(xx, yy):
    p, q = curv2rect(x1,x2,x3,x4,y1,y2,y3,y4,x,y)
    pp.append(p)
    qq.append(q)
pp = N.array(pp)
qq = N.array(qq)

# Plot
xp = [x1, x2, x3, x4, x1]
yp = [y1, y2, y3, y4, y1]
P.subplot(211)
levels = N.array([-10, 0, 1, 10.])
o = P.contourf(xx.reshape(nxy), yy.reshape(nxy), pp.reshape(nxy), levels=levels)
P.colorbar(o)
P.plot(xp, yp, 'k')
P.subplot(212)
o = P.contourf(xx.reshape(nxy), yy.reshape(nxy), qq.reshape(nxy), levels=levels)
P.colorbar(o)
P.plot(xp, yp, 'k')
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
P.savefig(figfile)
P.close()

