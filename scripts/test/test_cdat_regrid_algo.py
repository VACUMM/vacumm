"""Testing CDAT regridding algorithm"""
from vcmq import create_grid2d, meshbounds, set_grid, code_file_name, P, N, MV2, rc, add_grid
from collections import OrderedDict

configs = OrderedDict(
    libcf=['linear'], 
    esmf=[
        'linear',
        'patch',
        'conservative',
        ]
)

# Input curved grid
nxi = 5
nyi = 4
xxi, yyi = N.meshgrid(N.arange(nxi)+.25, N.arange(nyi)-.25)
for j in xrange(nyi):
    xxi[j,:] -= j*0.5
    #yyi[j,:] += j
for i in xrange(nxi):
    yyi[:,i] += i*0.5
gridi = create_grid2d(xxi,yyi) # input cdms grid
xxib,yyib = meshbounds(xxi,yyi) # coordinates of cell corners

# Output curved grid
nxo = 7
nyo = 7
xxo, yyo = N.meshgrid(N.arange(nxo)+.5, N.arange(nyo)-.5)
grido = create_grid2d(xxo, yyo) # output cdms grid
xxob,yyob = meshbounds(xxo,yyo) # coordinates of cell corners

# Input field
vari = MV2.array(N.arange(nyi*nxi).reshape(nyi,nxi))+10.
vari[1,1:3] = 100
vari[:2,1:3] = MV2.masked
set_grid(vari, gridi) # set grid and axes
#gridi.setMask(vari.mask)

# Define plot function
figfile = code_file_name(ext=False)+'_%(ifig)i.png'
#'%(tool)s_%(method)s.png'
figfiles = []
rc('font',size=9)
kw = dict(vmin=vari.min(),vmax=vari.max())
    
# Define logger
logfile = code_file_name(ext='log')
f = open(logfile, 'w')
def log(f, text): # logger
#    print text
    print >>f, text

# Loop on methods and tools
for tool, methods in configs.items():
    for method in methods:
        log(f, '-'*60)
        log(f,  ('%s/%s'%(tool,method)).upper())
        diag = {'dstAreaFractions': None,'dstAreas':
                None,'srcAreaFractions':None,'srcAreas':None} # output diags
        
        # Regridding
        varo = vari.regrid(grido, tool=tool, method=method, 
            diag=diag, coordSys='cart')
        
        # Ajustment for conservative methode
        frac = diag['dstAreaFractions']
        if method=='conservative': # divide by dstAreaFractions for conservative
            log(f, ' frac: %s'%frac)
            mask = frac==0.
            frac[mask] = 1.
            varo[:] /= frac
            varo[:] = MV2.masked_where(mask, varo, copy=0)
            log(f, ' dstareas: %s'%diag['dstAreas'])
            
        log(f, ' varo: %s'%varo)
        
        # Plot
        P.figure(figsize=(6,3))
        P.subplot(121).set_aspect(1)
        P.pcolormesh(xxib,yyib,vari.asma(),**kw)
        #P.colorbar()
        P.title('Original')
        kwg = dict(alpha=1, linewidth=1., centers=True)
        add_grid(gridi, color=(0,0,.2), marker='o', **kwg)
        add_grid(grido, color=(.2,0,0), marker='+', markersize=10,**kwg)
        P.axis('image')
        axis = P.axis()
        P.subplot(122).set_aspect(1)
        P.pcolormesh(xxob,yyob,varo.asma(),**kw)
        #P.colorbar()
        P.title('%(tool)s / %(method)s'%locals())
        add_grid(gridi, color=(0,0,.2), marker='o', **kwg)
        add_grid(grido, color=(.2,0,0), marker='+',markerlinewidth=1,markersize=8,**kwg)
        P.axis(axis)
        ifig = len(figfiles)
        ff = figfile%vars()
        P.tight_layout()
        P.savefig(ff)
        figfiles.append(ff)
        P.close()
f.close()

