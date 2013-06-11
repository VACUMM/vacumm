"""Compare CDAT regridding speed with rectangular and rectangular grids"""

config = {
   'esmf':['linear', 'patch', 'conserv'], 
   'libcf':['linear'], 
}

# Imports
from vcmq import *
from vacumm.misc.grid import rotate_grid
from time import time

# Input
nx = ny = 300
vari = MV2.array(N.arange(nx*ny*1.).reshape(ny, nx))
gridi = create_grid2d(vari.getAxis(1)[:]*50/nx,  vari.getAxis(0)[:]*50/nx)
set_grid(vari, gridi)

# Output grid
gridor = create_grid2d(vari.getAxis(1)[:]*0.09*50/nx,  
    vari.getAxis(0)[:]*0.09*50/nx)
gridoc = rotate_grid(gridi, 30)

# Log
logfile = code_base_name(ext='log')
if os.path.exists(logfile): os.remove(logfile)
f = open(logfile, 'w')
print >>f, 'NY=%(ny)i, NX=%(nx)i'%locals()

# Loop on methods
for tool, methods in config.items():
    for method in methods:
        print tool.upper(), method
        print >>f, tool.upper(), method
        
        print ' rect'
        t0 = time()
#        varo = vari.regrid(gridor, tool=tool, method=method)
        r = CDATRegridder(vari, gridor, tool=tool, method=method)
        t1 = time()
        dt = t1-t0
        varo = r(vari)
        print >>f,  ' rect: ini=%5.1f s + reg=%5.2f s %s'%(dt, time()-t1, psinfo())
        print ' curved'
        t0 = time()
        r = CDATRegridder(vari, gridoc, tool=tool, method=method)
        t1 = time()
        dt = t1-t0
        varo = r(vari)
        print >>f,  ' curv: ini=%5.1f s + reg=%5.2f s %s'%(dt, time()-t1, psinfo())
        
f.close()

