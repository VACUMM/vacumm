"""Test :func:`~vacumm.misc.atime.interp_clim`"""
from vcmq import (code_file_name, interp_clim, MV2, N, create_time, lindates, curve,
    strftime)

# Original clim
N.random.seed(0)
s = N.resize(N.sin(N.linspace(0, 1, 13)[:12]*2*N.pi), (2, 12)).T
clim = MV2.array(s, fill_value=1e20)
p = curve(clim[:, 0], 'o-', show=False, subplot=211, title='Original climatology',
    xmin=-.5, xmax=11.5, xticks=range(12),
    xticklabels=[strftime('%b', '2000-%i'%i) for i in range(1, 13)])

# Target times
times = lindates('2000-01-01', '2001-12-31', 5, 'day')

#  Interpolations
for i, method in enumerate(('linear', 'cubic', )):
    climo = interp_clim(clim, times, method=method)
    c = curve(climo[:, 0], 'o-', color='gr'[i], show=False, label=method.title(),
        subplot=212, title='Interpolated climatology', legend=True,
        markersize=2,
        tight_layout=True)
c.savefig(code_file_name())
c.show()

