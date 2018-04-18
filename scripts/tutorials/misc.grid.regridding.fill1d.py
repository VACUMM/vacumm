# -*- coding: utf8 -*-
import cdms2, MV2, pylab as P, curve2
from vcmq import data_sample, fill1d

# Lecture du niveau de la mer horaire
f = cdms2.open(data_sample('mars3d.t.nc'))
xe = f('xe')
f.close()
xe.long_name = 'Original'

# On crée des trous
# - petits
xe[:4] = MV2.masked
xe[12:16] = MV2.masked
# - gros
xe[40:46] = MV2.masked

# On rempli les petits trous (5 heures max) par interpolation cubique
xef = fill1d(xe, method='cubic', maxgap=5)
xef.long_name = 'Rempli'
xef[:] += xef.max()/5. # on décalle pour les plots

# Plots
P.rc('font', size=9)
curve2(xe, show=False, linewidth=1.5, figsize=(6, 3), top=.88, bottom=.15)
curve2(xef, show=False, linewidth=1.5, title='Niveau de la mer')
P.legend()
P.rcdefaults()
