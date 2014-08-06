"""Test the fortran function :f:func:`extrap1d`"""
from vcmq import N, P,meshcells, minmax
from vacumm.misc.grid._interp_ import extrap1d


mv = 999.
vari = N.zeros((5,5))+mv
vari[1,2:4] = [2,3]
vari[2,1] = 1
vari[2,3] = 3
vari[3,3:] = [3,4]
vari[4,:2] = [0,1]
vari = N.asfortranarray(vari)

varo0 = extrap1d(vari, mv, 0)
varop1 = extrap1d(vari, mv, 1)
varom1 = extrap1d(vari, mv, -1)
varop2 = extrap1d(vari, mv, 2)

result = [
    ('AssertTrue', N.allclose(vari, varo0)), 
    ('AssertTrue', N.allclose(varop1[:, -1], [999., 3., 3., 4., 1.])), 
    ('AssertTrue', N.allclose(varom1[:, 0], [999., 2., 1., 3., 0.])), 
    ('AssertTrue', N.allclose(varop1[:, -1], varop2[:, -1])), 
    ('AssertTrue', N.allclose(varom1[:, 0], varop2[:, 0])), 
    ]
    
