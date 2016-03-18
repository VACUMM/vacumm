"""Test :func:`~vacumm.misc.grid.misc.get_distances`"""

from vcmq import code_file_name, P, os, N
from vacumm.misc.grid import get_distances

result = []

# Input points
N.random.seed(0)
xa = N.random.uniform(-6., -4., (2, 3))
ya = N.random.uniform(47., 48., (2, 3))
xb = N.random.uniform(-5., -3., (2, 3))
yb = N.random.uniform(47.5, 48.5, (2, 3))

# Haversine
# - full
dhf = get_distances(xa, ya, xb, yb, mode='haversine')
result.append(('assertEqual', [dhf.shape, xa.shape+xb.shape]))
result.append(('assertEqual', [int(dhf.mean()/1e3), 97]))
# - pairwise
dhp = get_distances(xa, ya, xb, yb, mode='haversine', pairwise=True)
result.append(('assertEqual', [dhp.shape, xa.shape]))
result.append(('assertEqual', [int(dhp.mean()/1e3), 92]))

