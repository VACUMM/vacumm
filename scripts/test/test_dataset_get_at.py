"""Test :meth:`~vacumm.data.misc.dataset.OceanDataset.get_temp` with ``at`` keyword and :meth:`~vacumm.data.misc.dataset.OceanDataset.get_temp_u`"""


# Imports
from vcmq import DS, data_sample

# Inits
ncfile = "menor.nc"

# Read data
ds = DS(data_sample(ncfile), 'mars', logger_level='critical')
uu = ds.get_u3d()
uv1 = ds.get_u3d(at='v')
uv2 = ds.get_u3d_v()

# Checks
assert uu.mean() != uv1.mean()
assert uv1.mean() == uv2.mean()
assert uv1.id == "u3d_v"
assert uv2.id == "u3d_v"
