"""Test :meth:`~vacumm.data.misc.dataset.Dataset.get_var` on multiple files and time selections"""

# Imports
from vcmq import DS, data_sample

## Glob pattern
#ds = DS(data_sample("mars2d.xyt.*.nc"), 'mars', logger_level='critical',
#        time=("2008-08-15 06", "2008-08-15 09"))
#assert ds.get_ssh().shape == (4, 50, 50)
#assert ds.get_ssh(time=("2008-08-15 06", "2008-08-15 07", 'co'),
#                  squeeze=True).shape == (50, 50)

# With dates
ds = DS(data_sample("mars2d.xyt.%Y%m%d%H.nc"), 'mars', logger_level='debug',
        time=("2008-08-15 06", "2008-08-15 09"), patfreq=(8, 'hours'))
assert ds.get_ssh(
    time=("2008-08-15 06", "2008-08-15 07", 'cc')).shape == (2, 50, 50)
