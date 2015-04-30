"""Test :class:`~vacumm.misc.color.anamorph_cmap` """

# Imports
from vcmq import P, anamorph_cmap, code_file_name, plot_cmap

# Array-like
cmapa_old = P.get_cmap('jet')
cmapa_new = anamorph_cmap(cmapa_old, .9)

# Function-like
cmapf_old = P.get_cmap('rainbow')
cmapf_new = anamorph_cmap(cmapf_old, [.05, .2, .8, .95])

# Plot
fig = P.figure(figsize=(6, 2))
kw = dict(fig=fig, figsize=None, aspect=.15, show=False, close=False)
P.subplot(221)
plot_cmap(cmapa_old, **kw)
P.subplot(222)
plot_cmap(cmapa_new, **kw)
P.subplot(223)
plot_cmap(cmapf_old, **kw)
P.subplot(224)
plot_cmap(cmapf_new, **kw)
P.tight_layout()
P.savefig(code_file_name(ext='png'))
P.close()

