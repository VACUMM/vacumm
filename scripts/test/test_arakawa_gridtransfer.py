"""Test :meth:`vacumm.data.misc.arakawa.CGrid.interp`"""

from vcmq import MV2, N, create_grid, create_dep, set_grid, map2, \
    CGrid, minmax, curve2, add_grid, ArakawaGridTransfer, set_loc

# Initial variable
grid = create_grid(N.arange(-7, 0.), N.arange(43, 50.))
dep = create_dep([-5000, -3000, -2000, -1000, -500, -300, -200, -100.])
varc = {}
varc['t'] = MV2.reshape(N.arange(grid.size()*len(dep))*1., (len(dep), )+grid.shape)
set_grid(varc['t'], grid)
varc['t'].setAxis(0, dep)
set_loc(varc['t'],'t')

# Arakawa managers
ag = CGrid()
gt = ArakawaGridTransfer('c', 'a')

# Interpolations within same C grid
for p in 'u', 'v', 'f', 'w':
    varc[p] = ag.interp(varc['t'], 't', p, mode='extrap')
    
# Interpolation to A grid
vara = {}
for p in 't', 'u', 'v', 'f', 'w':
    vara[p] = gt.interp(varc[p], mode='extrap')

# For unittest
result = []
for p in 't', 'u', 'v', 'f', 'w':
    result.append(('assertTrue', N.ma.allclose(varc['t'], vara[p])))
