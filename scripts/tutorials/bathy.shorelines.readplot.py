from vcmq import map2, ocean, P, N, Histolitt
from _geoslib import Polygon


# %% Create a map manually (or pass ``m="auto"`` to ``coast.plot``)
m = map2(lon=(-5.2,  -3.4), lat=(47.6, 48.8), show=False, res=None,
         fillcontinents=False, drawcoastlines=False, figsize=(5.5, 4),
         bgcolor=ocean, left=.12, top=.9, proj='merc')

# %% Read Histolitt shape file
coast = Histolitt(clip=m.map)

# %% Plot filled polygons
coast.plot(color='0.8', linewidth=.5, edgecolor='.4',  m=m, show=False,
           title='Histolitt shoreline')

# %% Add some polygon areas
for poly, tpoly in zip(coast.shapes, coast.transform(m).shapes):
    area = tpoly.area() * 1e-6
    if area > 10:
        m.add_text(poly.boundary[:, 0].mean(), poly.boundary[:, 1].mean(),
                   "{:.0f}$km^2$".format(area), shadow=True, weight='bold',
                   transform='data', offset=(0, 10), ha='center', color='w',
                   bbox={'boxstyle': 'round', 'facecolor': 'tab:blue',
                         'alpha': .8, 'edgecolor': 'k'})
