# Fabrication des fausses donnees
import numpy as N
# - axe de temps
t = N.linspace(0., 100., 200)
# - observations
obs = N.sin(t)
# - modele
nm = 50
mod = N.resize(obs, (nm, len(t)))
mod += N.random.uniform(-1,  1, mod.shape)+.5

# Fabrication
from genutil.statistics import rms
bias = (mod-obs).mean(axis=1)
crms = rms(mod, N.resize(obs, mod.shape), centered=1, axis=1)
stdmod = mod.std(axis=1)
stdref = obs.std()

# Plot
from vacumm.misc.plot import target
target(bias, crms, stdmod, stdref, 
    savefigs=__file__,  savefigs_pdf=True, show=False)
