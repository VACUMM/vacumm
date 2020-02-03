"""Test :class:`~vacumm.misc.color.whiten` """

# Imports
from vcmq import P, whiten, plot_cmap, desaturate

cmap_name = "cmocean_thermal"
cmap_old = P.get_cmap(cmap_name)


# Max
cmap_10 = whiten(cmap_old, 1)

# Medium
cmap_05 = whiten(cmap_old, 0.5)

# Min
cmap_00 = whiten(cmap_old, 0)

# Plot
fig = P.figure(figsize=(6, 2))
kw = dict(fig=fig, figsize=None, aspect=.15, show=False, close=False)
P.subplot(221)
plot_cmap(cmap_old, title='Original', **kw)
P.subplot(222)
plot_cmap(cmap_00, title='f = 0.0', **kw)
P.subplot(223)
plot_cmap(cmap_05, title='f = 0.5', **kw)
P.subplot(224)
plot_cmap(cmap_10, title='f = 1.0', **kw)
P.tight_layout()
P.figtext(.5, 1, 'Whiten '+cmap_name, va='top', ha='center', size=12)
P.savefig(__file__[:-2]+'png')
P.close()

