import time,datetime,cdtime
from vacumm.misc.atime import is_datetime,is_cdtime

# Le temps de base : time
# - heure locale
mytime =  time.localtime()
print mytime 
#  -> (2007, 11, 12, 17, 5, 20, 0, 316, 0)
year = mytime[0]
# - chaine de caracteres
print time.asctime()
#  -> Mon Nov 12 17:08:13 2007

# Le temps manipulable (creation, transformation) : datetime
# !! ne fonctionne pas pour les dates avant 1900 !!
mytime = datetime.datetime(2000,10,1,2)
print mytime.day,mytime.second
#  -> 1 0
# - incrementer
mytime2 = mytime + datetime.timedelta(1,1) # (day,second)
print mytime2.day,mytime.second
#  -> 2 1

# Temps cdtime
#  ! Module de temps officiel d'vacumm !
#  Voir le tutoriel (*@\ref{lst:misc.time.bases.cdtime}@*) pour plus d'infos
ctime = cdtime.comptime(2000,10)
print mytime.year,mytime.month
#  -> 2000 10

# Verification des types
print is_datetime(mytime),is_cdtime(mytime)
# -> True False
print is_datetime(mytime),is_cdtime(ctime)
# -> True True
