"""Test the fortran function :f:func:`interp1dxx`"""
from vcmq import N, P,meshcells
from vacumm.fortran.interp import interp1dxx


nx = nyi = 10
mv = 1.e20
u, v = N.mgrid[-3:3:nx*1j, -3:3:10j]-2
vari = N.ma.asarray(u**2+v**2)
vari.set_fill_value(mv)
yi = N.linspace(-1000.,0., nyi)
yo = N.linspace(-1200, 100, 30.)
vari[int(nx/3):int(2*nx/3), int(nyi/3):int(2*nyi/3)] = N.ma.masked
x = N.arange(nx)
dyi = (yi[1]-yi[0])*0.49
dyo = (yo[1]-yo[0])*0.49
yyi = N.resize(yi, vari.shape)+N.random.uniform(-dyi, dyi, vari.shape)
yyo = N.resize(yo, (nx, len(yo)))+N.random.uniform(-dyo, dyo, (nx, len(yo)))
yyib, xxib  = meshcells(yyi, x)
yyob, xxob  = meshcells(yyo, x)

varon = N.ma.masked_values(interp1dxx(vari.filled(), yyi, yyo, mv, 0, extrap=0), mv)
varol = N.ma.masked_values(interp1dxx(vari.filled(), yyi, yyo, mv, 1, extrap=0), mv)
varoh = N.ma.masked_values(interp1dxx(vari.filled(), yyi, yyo, mv, 3, extrap=0), mv)

kw = dict(vmin=vari.min(), vmax=vari.max())
axlims = [x[0], x[-1], yo[0], yo[-1]]
P.figure(figsize=(8, 8))
P.subplot(221)
P.pcolor(xxib, yyib, vari)
P.axis(axlims)
P.title('Original')
P.subplot(222)
P.pcolor(xxob, yyob, varon, **kw)
P.axis(axlims)
P.title('Nearest1dxx')
P.subplot(223)
P.pcolor(xxob, yyob, varol, **kw)
P.axis(axlims)
P.title('Linear1dxx')
P.subplot(224)
P.pcolor(xxob, yyob, varoh, **kw)
P.axis(axlims)
P.title('Hermit1dxx')
P.tight_layout()
