"""Test function :func:`~vacumm.misc.grid.resol`"""

from vcmq import create_lon, create_lat, create_grid, rotate_grid, assert_allclose
from vacumm.misc.grid.misc import resol

lon1d = create_lon((0, 10., 1.))
lat1d = create_lat((43, 50, 1.))
rgrid = create_grid(lon1d, lat1d)
cgrid = rotate_grid(rgrid, 45)
lon2d = cgrid.getLongitude()
lat2d = cgrid.getLatitude()

assert_allclose(resol(lon1d, cache=False), 1.)
assert_allclose(resol(lon1d, meters=True, cache=False), 78626.2617245)

assert_allclose(resol(lon1d, lat=0, cache=False), 1.0)
assert_allclose(resol(lon1d, lat=0, meters=True, cache=False), 111195.031364)

assert_allclose(resol(rgrid, cache=False), (1.0, 1.0))
assert_allclose(resol(rgrid, meters=True, cache=False), (77242.051980165023, 111195.03136431401))

assert_allclose(resol((lon1d, lat1d), cache=False), (1.0, 1.0))
assert_allclose(resol((lon1d, lat1d), meters=True, cache=False), (77242.051980165023, 111195.03136431401))

assert_allclose(resol(cgrid, cache=False), (1.0, 1.0))
assert_allclose(resol(cgrid, meters=True, cache=False), (95933.976255,  95933.976255)) #(94804.663264739473, 94804.663264739473))


