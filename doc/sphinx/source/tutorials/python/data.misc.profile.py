#!/usr/bin/env python

import datetime
import cdms2, numpy, pylab
import vacumm.data.misc.profile as P

# Turn off netcdf4 file writting
cdms2.setNetcdfShuffleFlag(0); cdms2.setNetcdfDeflateFlag(0); cdms2.setNetcdfDeflateLevelFlag(0)

##########
# NOTES:
##########
# Profile:          Single profile object having one or more variables at a specific location and time
# Profiles:         Collection of Profile object with a save method suitable for ProfilesDataset
# ProfilesDataset:  Profiles dataset class with plot features
##########


##########
# 1) Load a profiles file using the generic Profiles class
##########
p = P.Profiles(
    # time, depth, latitude and longitude must be provided if their name are not those expected
    variables_map={
        'time':('TIME',),
        'depth':('DEPH',), # Force using the DEPH variable as depth
        'latitude':('LATITUDE',),
        'longitude':('LONGITUDE',),
        'temperature':('TEMP',),
    },
    # The variable argument specify which variables will be used when using the load/save methods.
    # Do not specify time, depth, latitude and longitude as they
    # are special variables which are always required
    variables=('temperature',),
    # We could also add a quality code filter:
    # qualities=(1,2)
)

from vacumm.config import data_sample
p.load(data_sample('data.misc.profile.nc'))

##########
# 2) Save profiles into a file suitable for profiles.py, that can be loaded by ProfilesDataset
##########
of = 'data.misc.profile.converted.nc'
p.save(of)
# We could also have specified the variables to write here:
# p.save(of, variables=('temperature',))

##########
# 3) Load the previously converted file using ProfilesDataset, which is more restrictive,
# it is expecting input file to have profile and level axis, time, latitude, longitude variables.
##########
p = P.ProfilesDataset(of)

# Plot the last two temperature profiles
pylab.figure()
p.plot_pro('temperature', -1, -2, plot_title='plot_pto')
pylab.savefig('data.misc.profile.plotpro.png')
pylab.figure()
p.plot_dist(
    plot_lon_min=-10, plot_lat_min=42,
    plot_lon_max=0, plot_lat_max=50,
    plot_title='plot_dist')
#from _geoslib import Polygon
#p.plot_dist(
#    polygons=(Polygon(numpy.array([[-10,42],[0,42],[0,50],[-10,50]])),),
#    plot_lon_min=-10, plot_lat_min=42,
#    plot_lon_max=0, plot_lat_max=50,
#    plot_title='plot_dist')
pylab.savefig('data.misc.profile.plotdist.png')
#pylab.show()


##########
# 4) Create a profile file from cdms variables
##########
# Populate a profile list with different number of levels for each profiles
npro = 5
profiles = []
for i in range(1,npro):
    profiles.append(
        P.Profile(
            platform_code='PRF%d'%i,
            datetime=datetime.datetime(2000, 01, 01, i),
            latitude=46+0.1*i,
            longitude=-8+0.1*i,
            depth=cdms2.createVariable(range(i))*-10.0,
            variables=dict(
                var1=cdms2.createVariable(numpy.random.rand(i)),
                var2=cdms2.createVariable(numpy.random.rand(i)))))
# Add another profile having only var1 variable
i += 1
profiles.append(
        P.Profile(
            platform_code='PRF%d'%i,
            datetime=datetime.datetime(2000, 01, 01, i),
            latitude=46+0.1*i,
            longitude=-8+0.1*i,
            depth=cdms2.createVariable(range(i))*-10.0,
            variables=dict(var1=cdms2.createVariable(numpy.random.rand(i)))))

# Note that variables must be specified, otherwise the default ('temperature','salinity') is used
# (feature used by merge_profiles.py, this must be fixed by a variables auto detection)
p = P.Profiles(profiles, variables=('var1', 'var2'))
p.save('data.misc.profile.fromscratch.nc')


