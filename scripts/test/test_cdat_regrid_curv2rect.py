"""Testing CDAT regridding from curvilinear grid to rectangular grid"""

config = {
    #'esmf':['conserv'],
   'esmf':[
        'linear',
        'patch',
        'conserv',
   ],
   'libcf':['linear'],
}
ncfile = 'swan.four.nc'
ncvar = 'HS'
nt = 10 ; nz = 30
nt = 1 ; nz = 1
xslice=slice(0, 100)
yslice=slice(0, 121)

# Imports
#from matplotlib import use ; use('Agg')
from vcmq import cdms2, bounds1d, bounds2d, data_sample, MV2, cp_props, N, minmax, code_file_name, P, os, psinfo, gc
#from cdms2.mvCdmsRegrid import CdmsRegrid
from time import time
from traceback import format_exc
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


#if rank==0: print 'read'
f = cdms2.open(data_sample(ncfile))
vari2d = f(ncvar)[0,yslice,xslice]
f.close()
gridi = vari2d.getGrid()
vari2d[:] += 3*vari2d.mean()
vari2d[:, -1] = MV2.masked
nyi,nxi = vari2d.shape

#if rank==0: print 'expand in time and depth'
vari = MV2.resize(vari2d, (nt, nz)+vari2d.shape)
cp_props(vari2d, vari)


#if rank==0: print 'grido'
loni = gridi.getLongitude()
lati = gridi.getLatitude()
xib, yib = bounds2d(loni, lati)
loni.setBounds(xib)
lati.setBounds(yib)
xi = loni.getValue()
yi = lati.getValue()
dx = N.diff(xi[0]).mean()
dy = N.diff(yi[:, 0]).mean()
xo = N.arange(xi.min()+10*dx, -30*dx+xi.max(), dx)
yo = N.arange(yi.min()-20*dy, yi.max()-20*dy, dy)
lono = cdms2.createAxis(xo)
lono.designateLongitude() ; lono.units= 'degrees_east'
lato = cdms2.createAxis(yo)
lato.designateLatitude() ; lato.units = 'degrees_north'
xob = bounds1d(lono) ; lono.setBounds(xob)
yob = bounds1d(lato) ; lato.setBounds(yob)
grido = cdms2.createRectGrid(lato, lono)
xmin, xmax = minmax(loni.asma(),lono)
ymin, ymax = minmax(lati.asma(), lato)
nyo,nxo = grido.shape
#print 'rank',rank
basefile = code_file_name(ext=False)
repfile = basefile+'.nt%(nt)s-nz%(nz)s-nyi%(nyi)s-nxi%(nxi)s.log'%locals()
if rank==0:
    if os.path.exists(repfile): os.remove(repfile)
    f = open(repfile, 'w')
    if size:
#        print 'MPI', size
        print >>f, 'MPI: NPROC=%i'%size
    print >>f, 'NT=%(nt)i, NZ=%(nz)i'%locals()
    print >>f, 'NYI=%(nyi)i, NXI=%(nxi)i'%locals()
    print >>f, 'NYO=%(nyo)i, NXO=%(nxo)i'%locals()

#if rank==0: print 'regridding...'
vmin = vari.min()
vmax = vari.max()
for tool, methods in config.items():
    for method in methods:
        if rank==0:
            t0 = time()
            print >>f, tool.upper(), method, ':'
            diag = {'dstAreaFractions': None}
        else:
            diag = None
        diag = {}
        try:
            varo = vari.regrid(grido, tool=tool, method=method, diag=diag)
            if rank==0: print >>f, ' ok'
        except:
            if rank==0:
                print >>f, ' failed'
                print >>f, format_exc()
            continue
        print 'diag:', diag
        if rank==0 and False:
            print >>f, tool.upper(), method, ':', '%5.1f'%(time()-t0), 'seconds', psinfo()
            frac = diag['dstAreaFractions']
            if frac is not None:
                mask = frac<=1.e-3
                frac[mask] = 1.
                frac = N.resize(frac, varo.shape)
                mask = N.resize(mask, varo.shape)
                varo[:] /= frac
                varo[:] = MV2.masked_where(mask, varo, copy=0)
#        del r
        gc.collect()
        if rank==0:
            print >>f, ' plot'
            P.figure(figsize=(12, 6))
            P.subplots_adjust(right=0.9)
            P.subplot(121)
            P.pcolormesh(xi, yi, vari[0,0].asma(),vmin=vmin,vmax=vmax)
            P.axis([xmin, xmax, ymin, ymax])
            P.colorbar()
            P.title('Original')
            P.subplot(122)
            P.pcolormesh(xo, yo, varo[0,0].asma(),vmin=vmin,vmax=vmax)
            P.axis([xmin, xmax, ymin, ymax])
            P.title(tool.upper()+' / '+method.upper())
            P.colorbar(extend='min')#cax=P.axes([0.92, 0.3, 0.02, 0.6]))
            figfile = basefile+'_%(tool)s_%(method)s.png'%vars()
            if os.path.exists(figfile): os.remove(figfile)
            P.savefig(figfile)
            P.close()
        del varo
if rank==0:print >>f, 'Done'
f.close()

