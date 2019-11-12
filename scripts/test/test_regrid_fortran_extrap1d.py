"""Test the fortran function :f:func:`extrap1d`"""
from utils import np
from vacumm.fortran.interp import extrap1d

vari = np.zeros((5,5)) + np.nan
vari[1, 2:4] = [2, 3]
vari[2, 1] = 1
vari[2, 3] = 3
vari[3, 3:] = [3, 4]
vari[4, :2] = [0, 1]

varo0 = extrap1d(vari, 0)
varop1 = extrap1d(vari, 1)
varom1 = extrap1d(vari, -1)
varop2 = extrap1d(vari, 2)

np.allclose(vari, varo0)
np.allclose(varop1[:, -1], [999., 3., 3., 4., 1.])
np.allclose(varom1[:, 0], [999., 2., 1., 3., 0.])
np.allclose(varop1[:, -1], varop2[:, -1])
np.allclose(varom1[:, 0], varop2[:, 0])
