from vacumm.misc.atime import strftime,strptime

# Lecture a partir d'une chaine de caracteres et d'un format
# (tappez strptime dans google)
mytime = strptime('1950-01-01 07:00:00','%Y-%m-%d %H:%M:%S')
# Verification partielle
print mytime.year,mytime.minute
# -> 1950 0

# On choisit la langue francaise 
import locale
locale.setlocale(locale.LC_ALL,'fr_FR')

# Ecriture sous un autre format
print strftime('%e %B %Y a %Hh%M',mytime)
# ->  1 janvier 1950 a 07h00
