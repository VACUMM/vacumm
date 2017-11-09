"""Test :meth:`~vacumm.misc.core_plot.Plot2D.plot_streamplot`"""

from vcmq import map2, data_sample, cdms2

f = cdms2.open(data_sample('mars2d.xyt.nc'))
u = f('u', slice(0, 1), squeeze=1)
v = f('v', slice(0, 1), squeeze=1)
f.close()

map2((u, v), fill=False, contour=False, streamplot=True, streamplot_density=3,
    streamplot_linewidth='modulus', streamplot_lwmod=3, streamplot_color='modulus',
    title='Tidal stream lines', colorbar_shrink=0.7,
    right=1)
