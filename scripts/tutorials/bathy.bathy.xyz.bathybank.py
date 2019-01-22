"""Manage random bathymetry datasets with :class:`~vacumm.bathy.bathy.XYZBathyBank`"""
from __future__ import print_function
import numpy as np,  os
from vcmq import XYZBathyBank

basename = __file__[:-2]

# %% Create fake bathymetries
# - generator
def gene_bathy(xc, yc, r, n):
    import numpy as np
    noise = np.random.random(n)
    a = np.random.random(n)*np.pi*2
    r = r*np.random.random(n)
    x = xc + r*np.cos(a)
    y = yc + r*np.sin(a)
    return np.asarray([x, y, np.exp(-(x-xc)**2/r**2-(y-yc)**2/r**2)*30.+noise]).transpose()
# - south
fsouth = basename + 'bathy_south.xyz'
np.savetxt(fsouth, gene_bathy(-5.1, 48.1, .15, 200))
# - north
fnorth = basename + 'bathy_north.xyz'
np.savetxt(fnorth, gene_bathy(-5.1, 48.3, .15, 100))
# - large
print(__file__)
flarge = basename + 'bathy_large.xyz'
print('flarge',flarge)
np.savetxt(flarge, gene_bathy(-5.2, 48., .4, 300))

# Put them in a bank
# - from scratch
bank_file = basename + 'bank.cfg'
if os.path.exists(bank_file):
    os.remove(bank_file)
bank = XYZBathyBank(bank_file)
# - add bathy files
bank.add(fsouth, id='south', long_name='South')
bank.add(fnorth, id='north', long_name='North')
bank.add(flarge) # id auto
# similar to:
# >>> bank += flarge
# - we check
print(bank)

# %% A few changes
# - id
bank['bathy.bathy.xyz.bathybank.bathy_large'].id = 'large'
bank['large'].update(long_name='Large')
# similar to:
# >>> bank.rename('bathy.bathy.xyz.bathybank.bathy_large', 'large')
# - transparency
bank['north'].transp = False
# - long name
bank['large'].long_name = 'Large'
bank['large']['transp'] = True

# %% Remove one
bank.remove('south')
#  similar to::
#  >>> bank -= 'south'
#  >>> bank -= bsouth
#  >>> del bank['south']

# %% Load the data
bsouth = bank['north'].load()

# %% Plot them all
bank.plot(size=15, map_proj='cyl', map_res=None, show=False)
