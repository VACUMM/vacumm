"""Test :func:`~vacumm.misc.grid.regridding.griddata"""

from vcmq import (griddata, MV2, code_file_name, N, create_grid, set_grid,
                  map2)

# Generate data
# - reference
xr = N.arange(20.) - 25
yr = N.arange(10.) + 43.
xxr, yyr = N.meshgrid(xr, yr)
zzr = (N.sin(xxr*N.pi/6)*N.sin(yyr*N.pi/6) + \
    N.exp(-((xxr-7.)**2+(yyr-7.)**2)/4.**2))*100.
zzr -= zzr.mean()
zzr = N.ma.asarray(zzr)
zzr[5:, 10:] = N.ma.masked
# - input at random locations
ij = N.unique((N.random.rand(150)*zzr.size).astype('i'))
xi, yi, zi = xxr.flat[ij], yyr.flat[ij], zzr.flat[ij]
zi = N.ma.resize(zi, (3, zi.size))
# - format
zi = MV2.array(zi, copy=False, id='sst')
taxis = zi.getAxis(0)
taxis.units = 'hours since 2000'
taxis.axis = 'T'
ggo = create_grid(xr, yr)
zzr = MV2.array(zzr)
set_grid(zzr, ggo)

# Call and plot
kw = dict(vmin=zzr.min(), vmax=zzr.max(),
          lon=(xr[0], xr[-1]), cmap_lum=0.7, linewidth=.5,
          lat=(yr[0], yr[-1]), show=False, colorbar=False)
m = map2(zzr, title='Original', subplot=221, figsize=(8, 5), **kw)
m.add_point(xi, yi, color='k')
for i, method in enumerate(['nearest', 'linear', 'cubic']):

    # Interp
    zo = griddata(xi, yi, zi, ggo, method=method)[0]

    # Plot
    m = map2(zo, title=method.title(), subplot=(2, 2, i+2), **kw)
    m.add_point(xi, yi, color='k')

m.savefig(code_file_name(ext='png'))
