# -*- coding: utf8 -*-
# On cherche directement 'Bre' pour Brest
from vacumm.tide.station_info import StationInfo
station = StationInfo('Bre')
# ->:
#Chargement de la station suivante :
#  Nom        : Ile de Brehat (Port-Clos)
#  Position   :  3.0°W / 48.9°N
#  Zone       : http://www.shom.fr/fr_page/fr_act_oceano/maree/zone6_9.htm
#  BM45       : 3.8
#  BM95       : 1.35
#  NM         : 5.89
#  PBM        : 0.1
#  PHM        : 11.6
#  PM45       : 8
#  PM95       : 10.4
#  ZERO_HYDRO : -5.4
#Definition des termes accessible avec .definitions()
# Loupé !
# On récupère la première station trouvée lors de l'initialisation
# Mais bon, on vérifie quand même
print station.attributes()
# -> ['igs', 'psmsl', 'uhslc', 'gloss', 'shom', 'legos', 'latitude', 
#       'longitude', 'zone', 'phm', 'pm95', 'pm45', 'nm', 'bm45', 
#        'bm95', 'pbm', 'zero_hydro']
print station.name, station.longitude
# -> Ile de Brehat (Port-Clos) -3.0

# On peut se servir de station pour continuer à chercher
# car le fichier est déjà chargé

# On affiche finalement toute les stations contenant 'bre'
print '-'*70
station.search('bre')
print '-'*70
# ->
#   Nom        : Ile de Brehat (Port-Clos)
#...
#   Nom        : Les Heaux-de-Brehat
 #...
#   Nom        : Brest

# Ok, on récupère Brest et uniquement Brest
# - en sélectionnant la station d'identifiant SHOM='Brest'
brest = station.find(shom='Brest', verbose=False)
print brest.name
# -> Brest
# - ou par sa position (station la plus proche)
brest = station.find(position=(4.5,48.4), verbose=False)
print brest.longitude, brest.latitude
# -> 2.36777777778 51.0347222222
