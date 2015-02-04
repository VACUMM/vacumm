"""Test the traditionnal CDAT regrid2 regridder"""

from vcmq import MV2, create_grid, meshbounds, P, add_grid, N, bounds1d, plot2d, savefigs,code_file_name
from regrid2 import Horizontal

# Input
nx, ny = 6, 4
vari = MV2.array(N.arange(nx*ny*1.).reshape(ny, nx), fill_value=1e20)
xi = vari.getAxis(-1)
xi[:] *= 2
yi = vari.getAxis(-2)
yi[:] *= 3
xi.designateLongitude()
yi.designateLatitude()
xi.setBounds(bounds1d(xi))
yi.setBounds(bounds1d(yi))
vari[1:2, 2:4] = MV2.masked
gridi = vari.getGrid()


# Output
grido = create_grid(xi[:]+2*2.5, yi[:]+3*1.5)
xo = grido.getLongitude()
yo = grido.getLatitude()
xo.setBounds(bounds1d(xo))
yo.setBounds(bounds1d(yo))
xxob, yyob = meshbounds(xo, yo)

# Regridding
varo, wo = vari.regrid(grido, tool='regrid2', returnTuple=1)

# Plot
kw = dict(fill='pcolor', contour=False, xhide=True, yhide=True,
    xticks=[], yticks=[], cmap='jet',
    colorbar=False, show=False)
P.figure(figsize=(6, 3))
p = plot2d(vari, subplot=131, title='Original', **kw)
add_grid(gridi)
add_grid(grido)
P.axis('image')
p = plot2d(varo, subplot=132, title='Regridded',  **kw)
add_grid(gridi)
add_grid(grido)
P.axis('image')
P.subplot(133)
P.pcolor(xxob, yyob, wo)
add_grid(grido)
P.title("Output weights")
P.tight_layout()
savefigs(code_file_name(),pdf=True, verbose=False)
P.close()



