"""Test :func:`~vacumm.diag.thermdyn.mixed_layer_depth` on MENOR"""

# Inits
ncfile = "menor.nc"

# Imports
from vcmq import DS, data_sample, mixed_layer_depth, rc, map2, os, code_file_name

# Read data
ds = DS(data_sample(ncfile),'mars', logger_level='critical')
temp = ds.get_temp(squeeze=True)
sal = ds.get_sal(squeeze=True)
depth = ds.get_depth(squeeze=True)
kz = ds.get_kz(squeeze=True)

# Compute MLD
kw = dict(depth=depth, format_axes=True)
mld = {}
mld['deltatemp'] = mixed_layer_depth(temp, mode='deltatemp',**kw)
mld['deltadens'] = mixed_layer_depth((temp,sal), mode='deltadens', **kw)
mld['kz'] = mixed_layer_depth(kz, mode='kz', **kw)
vmin = min([v.min() for v in mld.values()])
vmax = max([v.max() for v in mld.values()])

# Plot it
rc('font', size=8)
for i,(mode, var) in enumerate(mld.items()):
    m = map2(var, fill='pcolormesh', nmax=20, vmin=vmin, vmax=vmax,
        subplot=(len(mld),1,i+1), figsize=(4.1,8),
        contour_linewidths=0.7, cmap='vacumm_rnb2_hymex', hspace=0.25, bottom=0.08,
        title='%%(long_name)s: mode = "%s"'%mode, show=False)
figfile = code_file_name(ext='png')
if os.path.exists(figfile): os.remove(figfile)
m.savefig(figfile)

