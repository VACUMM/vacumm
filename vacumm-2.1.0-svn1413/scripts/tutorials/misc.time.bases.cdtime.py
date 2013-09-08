import cdtime

# Creer un objet 'comptime' : temps absolu
# - en specifiant tout (annee,mois,jour,heure,minute,seconde)
ctime = cdtime.comptime(2000,1,1,0,0,0)
# - ou seulement les premiers arguments
ctime = cdtime.comptime(2000)
# - a partir d'une chaines de caracteres
ctime = cdtime.s2c('2000-1-1 0')
# - on verifie
print ctime
#  -> 2000-1-1 0:0:0.0
# - on verifie des valeurs numeriques
print ctime.year,ctime.day
#  -> 2000 1

# Creer un objet 'reltime' : temps relatif
# - en specifiant explicitement (valeur, unites CF)
rtime = cdtime.reltime(50, 'years since 1950')
print '%s | %s | %s' %(rtime,rtime.value,rtime.units)
#  -> 50.000000 years since 1950 | 50.0 | years since 1950
# - a partir d'une chaines de caracteres et d'unites
rtime = cdtime.s2r('2000-1-1 0','years since 1950')
print rtime.value
#  -> 50.0

# Operations
# - soustraction/addition
print ctime.add(1,cdtime.Year),'|',rtime.add(-1,cdtime.Year)
#  -> 2001-1-1 0:0:0.0 | 49.00 years since 1950
# - conversions
rtime2 = ctime.torel('days since 2000')
ctime2 = rtime.tocomp().add(1,cdtime.Year)
# - comparaison
print rtime2 == rtime
#  -> True
print ctime2 <= ctime
#  -> False

# Verification des types
from vacumm.misc.atime import is_comptime,is_reltime
print is_comptime(ctime),is_reltime(rtime)
