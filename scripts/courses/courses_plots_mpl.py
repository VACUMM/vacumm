#!/usr/bin/env python
# -*- coding: utf8 -*- 
""" Plots using matplotlib

Examples from : http://matplotlib.org/gallery.html


For Matlab users, it can be very friendly ...

"""

# Imports

# WARNING: These two lines are necessary for displaying figures !
import matplotlib
matplotlib.use('qt4agg')

import matplotlib.pyplot as P
import matplotlib.cbook as cbook
import matplotlib.mlab as mlab
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

import numpy as N


# -----------------------------------------------------------------------------------------------------------
print 10*'-'+' ... plot et al. ... '+10*'-'

# - A figure
P. figure()

# -----------------------------------------------------------------------------------------------------------
# - plot a line
x = N.linspace(0, 10)
line, = P.plot(x, N.sin(x), '--', linewidth=2)

dashes = [10, 5, 100, 5] # 10 points on, 5 off, 100 on, 5 off
line.set_dashes(dashes)

P.show()
# => Practice: change figure properties - line color, xlim, ylim, ...

# -----------------------------------------------------------------------------------------------------------
# - plot a line with dates
datafile = cbook.get_sample_data('aapl.csv', asfileobj=False)
print ('loading %s' % datafile)
r = mlab.csv2rec(datafile)

r.sort()
r = r[-30:]  # get the last 30 days


# first we'll do it the default way, with gaps on weekends
fig, ax = P.subplots()
ax.plot(r.date, r.adj_close, 'o-')
fig.autofmt_xdate()
P.show()

# -----------------------------------------------------------------------------------------------------------
# - contour / pcolor


# => Practice: make these smaller to increase the resolution
dx, dy = 0.05, 0.05

# generate 2 2d grids for the x & y bounds
y, x = N.mgrid[slice(1, 5 + dy, dy),
                slice(1, 5 + dx, dx)]

z = N.sin(x) ** 10 + N.cos(10 + y * x) * N.cos(x)

# x and y are bounds, so z should be the value *inside* those bounds.
# Therefore, remove the last value from the z array.
z = z[:-1, :-1]
levels = MaxNLocator(nbins=15).tick_values(z.min(), z.max())

# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = P.get_cmap('PiYG')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

P. figure()
P.subplot(2, 1, 1)
# pcolormesh ... same as pcolor but faster
im = P.pcolormesh(x, y, z, cmap=cmap, norm=norm)
P.colorbar()
# set the limits of the plot to the limits of the data
P.axis([x.min(), x.max(), y.min(), y.max()])
P.title('pcolormesh with levels')

P.subplot(2, 1, 2)
# contours are *point* based plots, so convert our bound into point
# centers
P.contourf(x[:-1, :-1] + dx / 2.,
             y[:-1, :-1] + dy / 2., z, levels=levels,
             cmap=cmap)
P.colorbar()
P.title('contourf with levels')


P.show()

# -----------------------------------------------------------------------------------------------------------
# - quiver

Y, X = N.mgrid[-3:3:100j, -3:3:100j]
U = -1 - X**2 + Y
V = 1 + X - Y**2
speed = N.sqrt(U*U + V*V)

P.figure()
P.streamplot(X, Y, U, V, color=U, linewidth=2, cmap=P.cm.autumn)
P.colorbar()
# => Practice: Try to display these streamline through a vector field using quiver.

P.show()


# # -----------------------------------------------------------------------------------------------------------
# - scatter
# Load a numpy record array from yahoo csv data with fields date,
# open, close, volume, adj_close from the mpl-data/example directory.
# The record array stores python datetime.date as an object array in
# the date column
datafile = cbook.get_sample_data('goog.npy')
price_data = N.load(datafile).view(N.recarray)
price_data = price_data[-250:] # get the most recent 250 trading days

delta1 = N.diff(price_data.adj_close)/price_data.adj_close[:-1]

# Marker size in units of points^2
volume = (15 * price_data.volume[:-2] / price_data.volume[0])**2
close = 0.003 * price_data.close[:-2] / 0.003 * price_data.open[:-2]

fig, ax = P.subplots()
ax.scatter(delta1[:-1], delta1[1:], c=close, s=volume, alpha=0.5)

ax.set_xlabel(r'$\Delta_i$', fontsize=20)
ax.set_ylabel(r'$\Delta_{i+1}$', fontsize=20)
ax.set_title('Volume and percent change')

ax.grid(True)
P.tight_layout()

P.show()

# - Funny scatter polar plot ...
Nb = 150
r = 2 * N.random.rand(Nb)
theta = 2 * N.pi * N.random.rand(Nb)
area = 200 * r**2 * N.random.rand(Nb)
colors = theta

P. figure()
ax = P.subplot(111, polar=True)
c = P.scatter(theta, r, c=colors, s=area, cmap=P.cm.hsv)
c.set_alpha(0.75)
# => Practice: change colomapr

# - To display figure on the screen
P.show()