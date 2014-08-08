"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_uvgbt` on MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import data_sample, DS, N, MV2, shapiro2d, map2, code_file_name

# Read data
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')
u, v = ds.get_uvgbt()

# Mask huge values
mask = (N.ma.sqrt(u**2+v**2)>2.).filled(True)
for var in u, v:
    var[:] = MV2.masked_where(mask, var, copy=False)

# Shapiro filter
for var in u, v:
    var[:] = shapiro2d(var)

# Plot
map2((u[0], v[0]), title='Geostrophic velocity', 
    quiver_samp=3, xymasked=False, quiver_norm=3, contour=False, 
    fill=False, figsize=(7, 7), quiver_linewidth=.3, quiver_width=0.003, 
    quiver_scale=10, right=1, colorbar_shrink=0.8, bottom=0.05, show=False, 
    savefig=code_file_name(), close=True)

