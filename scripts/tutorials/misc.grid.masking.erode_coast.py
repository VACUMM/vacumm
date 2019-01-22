import MV2, numpy as N, matplotlib.pyplot as P
from vcmq import erode_coast, add_key

# %% Initial filled
xx, yy = N.indices((50, 100), 'f')
x0 = y0 = 25
dxy = 30
var = MV2.exp(-((xx-x0)**2+(yy-y0)**2)/dxy**2)

# %% Masks
# - reference
mask = var.filled()>.9
mask[:, 50:] = True
mask[15:35, 65:85] = False
# - variable
var[:] = MV2.masked_greater(var, .8)
var[:, 50:] = MV2.masked

# %% Erode
vare = erode_coast(var, mask)

# %% Plots
P.figure(figsize=(6, 9))
P.subplots_adjust(hspace=.2)
P.subplot(311)
P.pcolormesh(var.asma())
P.title('Original variable') ; add_key(1, color='w')
P.subplot(312)
P.pcolormesh(mask)
P.title('Reference mask') ; add_key(2, color='w')
P.subplot(313)
P.pcolormesh(vare.asma())
P.title('With coastal erosion') ; add_key(3, color='w')
