# %% False data
import numpy as N
# - time axis
t = N.linspace(0., 100., 200)
# - observations
obs = N.sin(t)
# - model
nm = 50
mod = N.resize(obs, (nm, len(t)))
mod += N.random.uniform(-3,  3, mod.shape)+.5

# %% Make stats
from genutil.statistics import rms
bias = (mod-obs).mean(axis=1)
crms = rms(mod, N.resize(obs, mod.shape), centered=1, axis=1)
stdmod = mod.std(axis=1)
stdref = obs.std()

# %% Plot
from vacumm.misc.plot import dtarget
dtarget(bias, crms, stdmod, stdref, colors='rms', show=False, units='m',
    scatter_alpha=0.5, sizes=40, savefigs=__file__, close=True)
