"""Generate a target diagram"""
import numpy as N, pylab as P
from genutil.statistics import rms
from vcmqm import dtarget

N.random.seed(0)

# Fake data
# - time axis
nt = 200
t = N.linspace(0., 100., nt)
# - observations
obs = N.sin(t)
# - model
nm = 50  # number of model states
mod = N.resize(obs, (nm, nt))
mod += N.random.normal(scale=0.2, size=mod.shape)
mod *= N.resize(N.random.uniform(.8, 1.6, size=nm), (nt, nm)).T
mod += .2

# Make stats
bias = (mod-obs).mean(axis=1)
crms = rms(mod, N.resize(obs, mod.shape), centered=1, axis=1).filled()
stdmod = mod.std(axis=1)
stdref = obs.std()

# Plots
P.figure(figsize=(7, 8))
P.subplot(211)
P.subplots_adjust(bottom=.05, hspace=.2)
P.title('Data')
for i in range(nm):
    P.plot(t, mod[i], color='tab:red', lw=.2, alpha=.2,
           label='Model states' if not i else '')
P.plot(t, obs, color='tab:blue', label='Observations')
P.legend()
P.subplot(212, anchor='S')
dtarget(bias, crms, stdmod, stdref, colors='rms', show=False, units='m',
        scatter_alpha=0.5, sizes=40, colorbar=True)
