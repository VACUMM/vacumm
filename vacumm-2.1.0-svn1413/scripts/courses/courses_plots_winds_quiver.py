#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Plot ARPEGE wind for selected period
"""

import matplotlib
matplotlib.use('qt4Agg')

import sys
import cdms2, cdtime
import MV2

import glob
import os
import numpy as N
import warnings
import matplotlib.pyplot as P

bcolorr='\033[91m'
bcolorg='\033[92m'
ecolor='\033[0m'

from vacumm.misc.grid import set_grid
from vacumm.misc.misc import MV2_concatenate 
from vacumm.misc.axes import create_lon, create_lat
from vacumm.misc.grid import curv_grid
from vacumm.misc.atime import ch_units, strftime
from vacumm.misc.io import ncread_files
from vacumm.misc.plot import savefigs, map2 as map
from vacumm.misc import auto_scale

# --------------------------------------------------------------------------


# -- Parameters to tune
year='2007'
rdeb = cdtime.comptime(2007,10,18)
rfin = cdtime.comptime(2007,10,25)
milo=-7
malo=-1
mila=44.5
mala=50


# -- Parameters 
rep_arpege='/home/coriolis_exp/spool/co01/co0123/co012302/co01230203/co0123020302_v3/best_estimate/'
rep_arpege = os.path.join(rep_arpege,year)
froot='METEOFRANCE_ARPEGE_%Y%m%dT00Z.nc'

# -- Read winds
u10 = ncread_files(os.path.join(rep_arpege,froot),'u10m',(rdeb,rfin))
v10 = ncread_files(os.path.join(rep_arpege,froot),'v10m',(rdeb,rfin))

print u10.shape


m=None
for it in N.arange(u10.shape[0]):# 

    P.figure()
    mod = MV2.sqrt(u10**2+v10**2)
    levels = auto_scale(mod,nmax=10, vmin=0., vmax=10.)
    kwarg=dict(m=m, lon=(milo,malo), lat=(mila,mala), nofill=True,  quiverkey_value=2, quiver_scale=35, 
                quiver_norm=3, contour=False, quiver_alpha=.9, 
                quiver_width=7., title='Wind - '+strftime('%d/%m/%Y - %H:%M:%S',cdtime.reltime(u10.getTime().getValue()[it],u10.getTime().units).tocomp()),   
                show=False, levels=levels, 
                colorbar_shrink=.7, right=1, quiver_samp=1, cmap='cmap_jete')
    # => Practice : change map parameters.
    
    print strftime('%d/%m/%Y - %H:%M:%S',cdtime.reltime(u10.getTime().getValue()[it],u10.getTime().units).tocomp())
    
    
    m=map((u10[it,0,:,:],v10[it,0,:,:]),**kwarg )
    # => Practice: instead of a quiver map, generate a map with current amplitude. 
    
    P.show()
    # savefigs('Wind_'+strftime('%Y%m%d_%H%M%S',cdtime.reltime(u10.getTime().getValue()[it],u10.getTime().units).tocomp())+'.png',dpi=200)
    
