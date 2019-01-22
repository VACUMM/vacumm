"""Test the :func:`~vacumm.misc.io.read_shapefile` function"""
from vcmqm import data_sample, read_shapefile, basic_proj

shpfile = data_sample("ne_110m_land/ne_110m_land")

# %% Basic
shapes, shape_type = read_shapefile(shpfile)
assert shape_type == 'polygon'
assert len(shapes) == 127

# %% Clip
shapes, shape_type = read_shapefile(shpfile, clip=[-10, 42, 10, 51.])
assert len(shapes) == 3

# %% Transform
shapes, shape_type = read_shapefile(shpfile, transform=basic_proj)
assert shapes[0].area() > 360 * 180
