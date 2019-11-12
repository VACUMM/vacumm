"""Test the fortran function :f:func:`cellerr1d`, :f:func:`cellerr1dx` et :f:func:`cellerr1dxx`"""
import numpy as np
from vacumm.fortran.interp import cellerr1d #, cellerr1dx, cellerr1dxx
# TODO: add test for cellerr1dx, cellerr1dxx

# %% Academic example pure 1d
# - data
yi = np.array([-2., 2.,  6. , 12., 15.])
vari = np.array(yi)
vari[-1] = np.nan
errm = np.array([1, 2, 1., 2., 3.])
yo = np.array([5.,  15])
yob = np.array([0., 10., 20.])
vari = vari.reshape(1, -1)
errm = errm.reshape(1, -1)
errl = np.ones(1)
# - interp
varo, erro = cellerr1d(vari, yi, yo, yob, errm, errl)
# - truth
errot = np.array([
    1/np.sqrt(
        1/(errm[0, 1]**2+np.abs(yo[0]-yi[1])*errl[0]) +
        1/(errm[0, 2]**2+np.abs(yo[0]-yi[2])*errl[0])),
    np.sqrt(errm[0, 3]**2+np.abs(yo[1]-yi[3])*errl[0])
    ])
varot = np.array([
    vari[0, 1]/(errm[0, 1]**2+np.abs(yo[0]-yi[1])*errl[0]) +
    vari[0, 2]/(errm[0, 2]**2+np.abs(yo[0]-yi[2])*errl[0]),
    vari[0, 3]/(errm[0, 3]**2+np.abs(yo[1]-yi[3])*errl[0])
    ])*errot**2
# - check
np.testing.assert_allclose(erro.ravel(), errot)
np.testing.assert_allclose(varo.ravel(), varot)

## 1dx
#yi = yi.reshape(1, -1)
#varox, errox = cellerr1d(vari, yi, yo, mv, errm, errl)
#np.testing.assert_allclose(erro, errox)
#np.testing.assert_allclose(varo, varox)
#
## 1dxx
#yo = yo.reshape(1, -1)
#varox, errox = cellerr1d(vari, yi, yo, mv, errm, errl)
#np.testing.assert_allclose(erro, errox)
#np.testing.assert_allclose(varo, varox)
