"""Test the :class:`~vacumm.misc.io.list_forecast_files` class"""
from vcmq import data_sample, list_forecast_files

# %% Local files
#assert len(list_forecast_files(data_sample('mars2d.xyt.2008081400.nc'))) == 0
#files = list_forecast_files(data_sample('mars2d.xyt.2008081500.nc'))
#assert len(files) == 1
#assert files[0] == data_sample('mars2d.xyt.2008081500.nc')
#assert len(list_forecast_files(data_sample('mars2d.xyt.20080815??.nc'))) == 3
print(list_forecast_files(data_sample('mars2d.xyt.%Y%m%d%H.nc'),
                               time=('2008-08-00', '2008-08-16')))
assert len(list_forecast_files(data_sample('mars2d.xyt.%Y%m%d%H.nc'),
                               time=('2008-08-00', '2008-08-16'))) == 3
#assert len(list_forecast_files(data_sample('mars2d.xyt.%Y%m%d%H.nc'),
#                               time=('2008-08-15 09', '2008-08-16'))) == 2
