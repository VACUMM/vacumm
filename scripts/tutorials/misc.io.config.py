# -*- coding: utf8 -*-
# On charge le fichier
from ConfigParser import SafeConfigParser
config = SafeConfigParser()
import os;print os.getcwd()
config.read('misc.io.config.in.ini')

# List des sections
print config.sections()

# On récupère les unités pour la SST
print config.get('sst', 'units')
#  -> m/s

# La latitude max de la section par défaut
print config.defaults()['lat_max']
#  -> 50.0
# Latitude max de la sst = celle par défaut
print config.getfloat('sst', 'lat_max')+1
#  -> 51.0

# Subsitutions
# - contenu substitué
print config.get('wind', 'name')
#  -> Wind on Iroise
# - contenu brut
print config.get('wind', 'name', raw=True)
# -> Wind on %(zone)s
# - contenu substitué par la force
print config.get('wind', 'name',  vars=dict(zone='for Britanny'))
#  -> Wind on Britanny



# On vire une section
config.remove_section('wind')
# On en crée une autre
print config.has_section('sealevel')
#  -> False
config.add_section('sealevel')
config.set('sealevel','name', 'Sea level')
config.set('sealevel','units', 'm')
print config.has_option('sealevel', 'name')
#  -> True

# On sauvegarde
fc = open('misc.io.config.out.ini', 'w')
config.write(fc)
fc.close()
