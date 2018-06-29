"""Test the fortran function :f:func:`cellerr1d`, :f:func:`cellerr1dx` and :f:func:`cellerr1dxx`"""
import numpy as N
from vacumm.misc.grid._interp_ import cellerr1d, cellerr1dx, cellerr1dxx

# Academic example pure 1d
# - data
yi = N.array([-2., 2.,  6., 12., 15.])
vari = N.ma.array(yi)
vari[-1] = N.ma.masked
errm = N.ma.array([1, 2, 1., 2., 3.])
yo = N.array([5.,  15])
mv = 1e20
vari = vari.reshape(1, -1).filled(mv)
errm = errm.reshape(1, -1).filled(mv)
errl = N.ones(1)
# - interp
varo, erro = cellerr1d(vari, yi, yo, mv, errm, errl)
# - truth
errot = N.array([
    1/N.sqrt(
        1/(errm[0, 1]**2+N.abs(yo[0]-yi[1])*errl[0]) +
        1/(errm[0, 2]**2+N.abs(yo[0]-yi[2])*errl[0])),
    N.sqrt(errm[0, 3]**2+N.abs(yo[1]-yi[3])*errl[0])
    ])
varot = N.array([
    vari[0, 1]/(errm[0, 1]**2+N.abs(yo[0]-yi[1])*errl[0]) +
    vari[0, 2]/(errm[0, 2]**2+N.abs(yo[0]-yi[2])*errl[0]),
    vari[0, 3]/(errm[0, 3]**2+N.abs(yo[1]-yi[3])*errl[0])
    ])*errot**2
# - check
N.testing.assert_allclose(erro.ravel(), errot),
N.testing.assert_allclose(varo.ravel(), varot)

# 1dx
yi = yi.reshape(1, -1)
varox, errox = cellerr1d(vari, yi, yo, mv, errm, errl)
N.testing.assert_allclose(erro, errox)
N.testing.assert_allclose(varo, varox)

# 1dxx
yo = yo.reshape(1, -1)
varox, errox = cellerr1d(vari, yi, yo, mv, errm, errl)
N.testing.assert_allclose(erro, errox)
N.testing.assert_allclose(varo, varox)
