# -*- coding: utf8 -*-
import cdms2, pylab as P
from vcmq import data_sample, Marigraph

# Lecture du niveau de la mer
f = cdms2.open(data_sample('tide.sealevel.BREST.mars.nc'))
sea_level = f('sea_level', time=('2006-10', '2006-10-07'))
f.close()

# Initialisation de l'objet marégraphique
mg = Marigraph(sea_level, verbose=True)

# Extraction du signal de maree
tide = mg.tide(tide_filter='demerliac')

# On peut aussi specifier le filtre avec :
mg.set_tide_filter('demerliac') # On peut aussi choisir 'godin'

# Calcul surcotes/decotes
cotes = mg.cotes()

# Calcul des pleines et basses mers
ref = 'mean'
highs = mg.highs(ref=ref)
lows = mg.lows(ref=ref)
zeros = mg.zeros(ref=ref)

# On plot
P.rc('font', size=10)
P.figure(figsize=(5, 5))
kwplot = dict(date_fmt='%d/%m', date_locator='day', show=False, hspace=.2,  left=.15)
# - tout sauf le signal d'origine
P.subplot(211)
mg.plot(orig=False, tide_color='k', title='Extremas et zeros', xhide=True, **kwplot)
# - seules les surcotes decotes
P.subplot(212)
mg.plot('cotes', title='Surcotes et decotes', **kwplot)
