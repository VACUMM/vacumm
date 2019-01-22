"""Merging of random bathymetry measurements"""

# %% Make fake bathymetries
from vcmq import N, P, XYZBathy,  XYZBathyMerger
# - generator
def gene_bathy(xc, yc, xr, yr, n=500, amp=30.):
    import numpy as np
    noise = np.random.random(n)
    a = np.random.random(n)*np.pi*2
    r = np.random.random(n)
    x = xc + xr*r*np.cos(a)
    y = yc + yr*r*np.sin(a)
    return np.asarray([x, y, np.exp(-(x-xc)**2/xr**2-(y-yc)**2/yr**2)*amp+noise])
# - creations
xyz1 = XYZBathy(gene_bathy(-5, 48.3, .3, .15))  # top right
xyz2 = XYZBathy(gene_bathy(-5.45, 48.3, .2, .1,  amp=10))  # top left
xyz3 = XYZBathy(gene_bathy(-5.45, 48.1, .45, .25,
                           n=1000, amp=20), transp=False)  # bot left
xyz4 = XYZBathy(gene_bathy(-5., 48.1, .14, .08, n=300))  # bot right
fxyz5 = __file__[:-2]+'xyz5.xyz'
N.savetxt(fxyz5, gene_bathy(-5.15, 48.2, .2, .1, amp=15).transpose())  # center

# %% Direct merging
xyz = xyz1 + xyz2 + xyz3 + xyz4 + fxyz5

# %% Plot
P.figure(figsize=(5., 8.5))
P.subplot(311)
P.rcParams['font.size'] = 9
P.subplots_adjust(bottom=.03, top=.97, hspace=.25)
kwplot = dict(show=False, colorbar=False, map_res=None,
              margin=0., map_autoresize=0,
              xmin=xyz.xmin, xmax=xyz.xmax, ymin=xyz.ymin, ymax=xyz.ymax)
xyz.plot(title='Fusion directe', **kwplot)

# %% Use a merger
# - init
merger = XYZBathyMerger()
# - direct add
merger += xyz1
merger.append(xyz2)
merger += xyz3
merger += xyz4
# - add from a file
merger += fxyz5
# - remove one dataset
merger -= xyz2
# - remove the last one
del merger[-1]
# - change some attributes
for i, xyz in enumerate(merger):
    xyz.long_name = 'Niveau : {}'.format(i)  # For the legends
    xyz.set_transp(False)  # Opacity

# %% Plot the merger
# - values
P.subplot(312)
merger.plot(mode='value', title='Merger: values', **kwplot)
# - cluster
P.subplot(313)
merger.plot(mode='cluster', size=10, title='Merger: cluster', marker='o',
            legend_loc='upper left', **kwplot)

# %% Get it as x, y z
merged_xyz = merger.get_xyz(long_name='Merged')

# %% Save it
merged_xyz.save(__file__[:-2]+'merged.xyz')
