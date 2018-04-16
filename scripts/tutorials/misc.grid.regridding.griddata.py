# On construit un grille reguliere
import numpy as N
from vacumm.misc.grid import meshbounds
xr = N.arange(20.)
yr = N.arange(10.)
xxr, yyr = N.meshgrid(xr, yr)
xrb, yrb = meshbounds(xr, yr)
zzr = (N.sin(xxr*N.pi/6)*N.sin(yyr*N.pi/6) + \
    N.exp(-((xxr-7.)**2+(yyr-7.)**2)/4.**2))*100.
zzr -= zzr.min()
vminmax=dict(vmin=zzr.min(), vmax=zzr.max())

# On construit un echantillon irregulier a partir du regulier
ij = N.unique((N.random.rand(50)*zzr.size).astype('i'))
xi, yi, zi = xxr.flat[ij], yyr.flat[ij], zzr.flat[ij]

# Interpolation sur la grille reguliere
# - natgrid
from vacumm.misc.grid.regridding import griddata
zirn = griddata(xi, yi, zi, (xr, yr), method='nat', ext=True, sub=6)
# - krigeage
zirk = griddata(xi, yi, zi, (xr, yr), method='carg')

# Plot de verif
import pylab as P
from vacumm.misc.plot import savefigs
P.figure(1, figsize=(6, 8))
P.subplots_adjust(hspace=.3, bottom=.05, top=.95, left=.06)
# - regulier
P.subplot(311)
P.pcolor(xrb, yrb, zzr, **vminmax)
P.xlim(xrb.min(), xrb.max()) ; P.ylim(yrb.min(), yrb.max())
P.title('Original')
stdref = zzr.std()
# - irregulier via natgrid
P.subplot(312)
P.pcolor(xrb, yrb, zirn, **vminmax)
P.plot(xi, yi, 'ko')
P.xlim(xrb.min(), xrb.max()) ; P.ylim(yrb.min(), yrb.max())
P.title('Natgrid (err = %02i%%)'%((zzr-zirn).std()*100/stdref))
# - irregulier via cargen
P.subplot(313)
P.pcolor(xrb, yrb, zirk, **vminmax)
P.plot(xi, yi, 'ko')
P.xlim(xrb.min(), xrb.max()) ; P.ylim(yrb.min(), yrb.max())
P.title('Minicargen (err = %02i%%)'%((zzr-zirk).std()*100/stdref))
savefigs(__file__)
P.close()
