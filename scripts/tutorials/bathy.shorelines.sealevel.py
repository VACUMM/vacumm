"""Interpolate mean sea level onto the shoreline"""
from vacumm.bathy.shorelines import GSHHS

# %% Britanny zone
zone = (-6.5, 47.2, -2, 49.05)

# %% Get the shoreline
tc = GSHHS(input='l', clip=zone)

# %% Interpolate
xyz = tc.bathy()

# %% Plot
xyz.plot(size=10, show=False)
