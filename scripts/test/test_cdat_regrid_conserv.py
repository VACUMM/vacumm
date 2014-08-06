"""Test CDAT conservative regridding"""
from vcmq import N, meshbounds, bounds1d, cdms2, MV2, rc, P, add_grid, rcdefaults, \
    create_lon, create_lat, savefigs, code_file_name
    
# Input grid
x0, y0, nx, ny, dx, dy = 0., 0., 20, 15, 5., 5.
x = x0+N.arange(0, nx*dx, dx)
y = y0+N.arange(0, ny*dy, dy)
xxbi, yybi = meshbounds(x, y)
xb = bounds1d(x)
yb = bounds1d(y)
lon = create_lon(x)
lat = create_lat(y)
lon.setBounds(xb)
lat.setBounds(yb)
gridi = cdms2.createRectGrid(lat, lon)

# Input data
vari = MV2.ones(gridi.shape)*50
vari[5:10, 7:13] = -10
#vari[:, -1] = MV2.masked # <<< THIS MAKES IT WORK!

vari.setAxisList([gridi.getLatitude(), gridi.getLongitude()])
vari.setGrid(gridi)

# Output grid
grido, xxbo, yybo = gridi, xxbi, yybi

# Regrid
diag = {'dstAreaFractions':None}
varo = vari.regrid(grido, tool='esmf', method='conservative', 
    diag=diag, coordSys='cart')
    
# Norm
frac = diag['dstAreaFractions']
mask = frac==0.
frac[mask]=1.
varo[:] /= frac
varo[:] = MV2.masked_where(mask, varo, copy=0)

# Plot
rc('font', size=9)
kw = dict(vmin=vari.min(), vmax=vari.max())
axis = [xxbi.min(), xxbi.max(), yybo.min(), yybo.max()]
P.figure(figsize=(7, 3))
P.subplot(121, aspect=1)
P.pcolormesh(xxbi, yybi, vari, **kw)
P.colorbar(shrink=0.7)
add_grid(gridi, color='0.5')
add_grid(grido)
P.axis(axis)
P.title('Original: max=%g min=%g'%(vari.max(), vari.min()))
P.subplot(122, aspect=1)
P.pcolormesh(xxbo, yybo, varo, **kw)
P.colorbar(shrink=0.7)
add_grid(gridi)
add_grid(grido, color='0.5')
P.axis(axis)
P.title('Regridded: max=%g min=%g'%(varo.max(), varo.min()))
P.tight_layout()
savefigs(code_file_name(), verbose=False)
P.close()
rcdefaults()

