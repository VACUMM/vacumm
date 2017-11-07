#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot the MPI domain decomposition"""

# source: /work/mdussauz/MARS/RUN_MARS/VILAGR/VILAGR-r1023/rank_1
# showmpi.py -o mpi.png -s -m /home1/caparmor/sraynaud/projects/VACUMM3/scripts/actimar/divers/bathy_verif.nc /home1/caparmor/sraynaud/projects/VACUMM3/scripts/actimar/divers/mpi.txt
# -o mpi.png ../data/mpi.nc ../data/mpi.txt

import argparse
import os

# Arguments
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("marsfile", help="MARS file")
parser.add_argument("mpifile", help="MPI domain decomposition file [default: %(default)s]", default="mpi.txt")
parser.add_argument("-q", "--quiet", help="do print information", action='store_true')
parser.add_argument("-t", "--title", help="Plot title", default="MPI domain decomposition")
parser.add_argument("-a", "--alpha", help="apha transparency [0->1]", default=0.3, type=float)
parser.add_argument("-o", "--outfig", help="Output figure file")
parser.add_argument("-d", "--nodisp", action='store_true', help="Do not show the figure")
options = parser.parse_args()

# Imports
import numpy as N, MV2
import matplotlib.pyplot as P
from matplotlib.ticker import FormatStrFormatter, FuncFormatter
from warnings import warn
from vacumm.misc import scalebox
from vacumm.bathy.bathy import plot_bathy
from vacumm.misc.grid import meshbounds, meshgrid
from vacumm.misc.color import cmap_srs
from vacumm.data import setup_dataset

# Read the MARS file
if not os.path.exists(options.marsfile):
    parser.error('Invalid MARS file name: '+options.marsfile)
ds = setup_dataset('mars', options.marsfile, log_level='error' if options.quiet else 'info')
bathy = ds.get_bathy().asma()
bathy[:] = N.ma.masked_less(bathy, -30., copy=False)
mask = N.ma.getmaskarray(bathy) #; del bathy
valid = ~mask
mask = mask.astype('i')
ny, nx = mask.shape

# Read the MPI file
if not os.path.exists(options.mpifile):
    parser.error('Invalid MPI file name: '+options.mpifile)
f = open(options.mpifile)
nb = int(f.readline().split()[0]) # Number of procs
nbx = int(f.readline().split()[0])
nby = int(f.readline().split()[0])
sh = (nbx, nby)
mpi = N.recfromtxt(f, skip_header=2, names=True, converters={6:lambda v: 1 if v=='T' else 0}, dtype='i')
f.close()

# Find inactive ocean points
bads = {}
for ib in xrange(len(mpi)):
    sel = (slice(mpi.jmin[ib], mpi.jmax[ib]+1), slice(mpi.imin[ib], mpi.imax[ib]+1))
    if not mpi.actif[ib] and valid[sel].any():
        bads[ib] = valid[sel].sum()
        submask = mask[sel]
        mask[sel] = N.where(submask==0, 2, submask) ; del submask
        mpi.actif[ib] = 2
mpi2d = mpi.reshape(nby, nbx)
badblocks = mpi2d.actif==2

# Infos
if not options.quiet:
    nbb = badblocks.sum()
    if nbb==0:
        print 'Found no inactive ocean points'
    else:
        print 'Found a total of %i inactive ocean point(s) in %i MPI block(s): '%((mask==2).sum(), nbb)

# Indices of MPI blocks
iib = N.empty((nby+1, nbx+1))
iib[:nby, :nbx] = mpi2d.imin[:, :nbx]-.5
iib[:nby, -1] = mpi2d.imax[:, -1]+0.5
iib[-1] = iib[-2]
jjb = N.empty((nby+1, nbx+1))
jjb[:nby, :nbx] = mpi2d.jmin[:nby]-.5
jjb[-1, :nbx] = mpi2d.jmax[-1]+0.5
jjb[:, -1] = jjb[:, -2]
jj, ii = N.indices((nby, nbx))

# Plot model mask
P.imshow(mask, origin='lower', interpolation='none', aspect=1,
    cmap='jet', vmin=0, vmax=2)

# Plot MPI blocks
P.pcolor(iib, jjb, N.clip(mpi2d.actif, 0, 1),
    shading='faceteded', edgecolors='k', linewidths=.8,
    alpha=options.alpha, cmap='gray', vmin=0, vmax=1)
if badblocks.any(): # Highlight suspect blocks
    for i, j in zip(ii[badblocks], jj[badblocks]):
        P.plot([iib[i, j], iib[i+1, j], iib[i+1, j+1], iib[i, j+1], iib[i, j]],
            [jjb[i, j], jjb[i+1, j], jjb[i+1, j+1], jjb[i, j+1], jjb[i, j]],
            '-', color=(1, 0, 0), lw=2)
P.grid('off')
P.xlabel('Grid indices along X')
P.ylabel('Grid indices along Y')
P.axis('image')
fmt = FuncFormatter(lambda v, p=None: '%i'%round(v))
P.gca().xaxis.set_major_formatter(fmt)
P.gca().yaxis.set_major_formatter(fmt)
P.title(options.title)
P.tight_layout()
if options.outfig:
    P.savefig(options.outfig)
    print 'Plot saved to: '+options.outfig
if not options.nodisp:
    P.show()

