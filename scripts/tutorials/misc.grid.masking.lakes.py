# Inits
from __future__ import print_function
from vcmq import Lakes, meshbounds, cmap_linear, cmap_srs, land
import numpy as N, pylab as P

P.figure(figsize=(5.5, 5))
P.subplots_adjust(bottom=.02, top=.92, left=.01, right=.98, hspace=.3)
cmap_mask = cmap_linear(((.6, .8, 1), land))
cmap_lakes = cmap_srs(('w', 'r', 'g', 'b'), stretch=0)

# Create a mask with lakes
mask = N.ones((20, 30), '?') # Land
mask[10:, 10:] = False  # Ocean
mask[10:15, 13:14] = True  # Ocean
mask[0, 0] = False      # Little lake
mask[5:15, 2:4] = False # Big lake
xxb, yyb = meshbounds(N.arange(mask.shape[1]*1.), N.arange(mask.shape[0]*1.))
xlim = (xxb.min(), xxb.max()) ; ylim = (yyb.min(), yyb.max())
P.subplot(221)
P.pcolor(xxb, yyb, mask.astype('i'), cmap=cmap_mask)
P.xlim(xlim) ; P.ylim(ylim) ; P.xticks([])  ;P.yticks([])
P.title('Initial mask')

# Find lakes
lakes = Lakes(mask)
m = lakes[0]
P.subplot(222)
P.pcolor(xxb, yyb, lakes.lakes(), cmap=cmap_lakes)
P.xlim(xlim) ; P.ylim(ylim) ; P.xticks([])  ;P.yticks([])
P.title('All lakes')

# The second lake
P.subplot(224)
P.pcolor(xxb, yyb, lakes.lakes(2), cmap=cmap_lakes, vmax=lakes.nlakes)
P.xlim(xlim) ; P.ylim(ylim) ; P.xticks([])  ;P.yticks([])
P.title('Second lake only')

# Ocean
P.subplot(223)
ocean = lakes.ocean()
P.pcolor(xxb, yyb, ocean, cmap=cmap_mask)
P.xlim(xlim) ; P.ylim(ylim) ; P.xticks([])  ;P.yticks([])
P.title('Mask ocean')
print('Success?', (ocean==lakes[0]).all())
