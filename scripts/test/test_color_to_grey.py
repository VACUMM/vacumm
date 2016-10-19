"""Test :class:`~vacumm.misc.color.to_grey` """

# Imports
from vcmq import P, to_grey, code_file_name, plot_cmap


cmap_name = "CMRmap"
cmap_old = P.get_cmap(cmap_name)


# Max
cmap_10 = to_grey(cmap_old, 1)

# Medium
cmap_05 = to_grey(cmap_old, 0.5)

# Min
cmap_00 = to_grey(cmap_old, 0)

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
P.savefig(code_file_name(ext='png'))
P.show()
P.close()

