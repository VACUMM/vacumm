"""Fill masked data on 1D arrays"""
from vcmqm import plt, MV2, cdms2, data_sample, fill1d, curve2

# Read hourly sea level
f = cdms2.open(data_sample('mars3d.t.nc'))
xe = f('xe')
f.close()
xe.long_name = 'Original'

# Create fake holes
# - small
xem = xe.clone()
xem[:4] = MV2.masked
xem[12:16] = MV2.masked
# - big
xem[40:46] = MV2.masked
xem.long_name = 'Masked'

# Fill small holes using cubic interpolation
xef = fill1d(xem, method='cubic', maxgap=5)
xef.long_name = 'Filled [max gaps = 5]'

# Plot
plt.rc('font', size=9)
curve2(xef, show=False, linewidth=4, color='0.8', figsize=(6, 3),
       top=.88, bottom=.15)
curve2(xe, 'o', show=False, color='k', markersize=3)
curve2(xem, show=False, linewidth=1, color='tab:red', title='Sea level')
plt.legend(loc='upper right')
plt.rcdefaults()
