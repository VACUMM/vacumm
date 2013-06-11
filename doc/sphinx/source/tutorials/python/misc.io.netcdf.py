# Ouverture
import cdms2
from vacumm.config import data_sample
f = cdms2.open(data_sample('mars2d.xyt.nc'))

# Lister les variables
print f.variables.keys()
#  -> ['bounds_lon', 'h0', 'v', 'xe', 'u', 'bounds_lat']

# Avoir des informations sur une variables sans la lire, via []
nt = f['xe'].shape[0]
print f['xe'].getTime().asComponentTime()[0:nt:nt-1]
#  -> [2008-8-15 0:0:0.0, 2008-8-15 23:0:0.0]

# Lire une selection de la variable
import cdtime
xe = f('xe', ('2008-8-15',cdtime.comptime(2008,8,15,12), 'cc'), 
    lon=slice(5,6), lat=(48.1, 48.5), squeeze=1)
print xe.shape
# -> (13, 29)
# squeeze a supprime l'axes des longitudes de dim 1

# Fermer le fichier lu
f.close()

# Definir la compression netcdf4
cdms2.setNetcdfShuffleFlag(1)
cdms2.setNetcdfDeflateFlag(1)
cdms2.setNetcdfDeflateLevelFlag(3)

# Creer un nouveau fichier
ncfile = 'misc-io-netcdf.nc'
import os
if os.path.exists(ncfile): os.remove(ncfile) 
f = cdms2.open('misc-io-netcdf.nc', 'w') # ouverture en ecriture
f.write(xe) # ecriture d'une variable
f.history = 'Created with '+__file__.encode('utf8') # attribut global
f.close() # fermeture
