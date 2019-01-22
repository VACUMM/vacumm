# -*- coding: utf8 -*-
"""Sigma coordinates conversions using :class:`~vacumm.misc.sigma.NcSigma`"""
from __future__ import print_function
import cdms2, os
import numpy as np
from vcmq import (data_sample, NcSigma, create_depth, netcdf4, regrid1d,
                  section2)

# %% Output interpolation depths
depths = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
    22,24,26,28,30,32,34,36,38,40,45,50,55,60,65,70,75,80,85,90,95,100,120,140,160])
depths = -depths[::-1] # croissantes et négatives

# %% Read data
f = cdms2.open(data_sample('mars3d.tsigyx.nc'))
data_in = f('TEMP') # T-YX (- = level)

# %% Compute input sigma depths
# - guess sigma type and class from file meta-data
sigma_class = NcSigma.factory(f)
# - initialise the converter
sigma_converter = sigma_class(copyaxes=True)
# - read what is needed and make conversion
depths_in = sigma_converter().filled()
f.close()

# %% Create the Z axis
depth_out = create_depth(depths)

# %% Interpolate to constant depths
data_out = regrid1d(data_in, depth_out, axi=depths_in, axis=1,
                    method='linear', extrap=1)

# %% Plot
kw = dict(show=False, vmin=10, vmax=14, xhide='auto', add_grid=True, ymax=0,
          cmap='cmocean_thermal')
section2(data_in[0, :, 10], yaxis=depths_in[0, :, 10], subplot=211,
         title='Sigma', **kw)
s = section2(data_out[0, :, 10], subplot=212, title='Z', **kw)

# Sauvegarde
outfile =__file__ [:-2] + '.nc'
if os.path.isfile(outfile):
  os.remove(outfile)
f2 = cdms2.open(outfile,'w')
f2.write(data_out)
f2.close()
print('Saved to', outfile)

