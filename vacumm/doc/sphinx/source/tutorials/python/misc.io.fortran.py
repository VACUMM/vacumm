from Scientific.IO.FortranFormat import FortranFormat,FortranLine

# Fichier ascii exemple
f = open('misc.io.fortran.dat', 'w')
f.write('   59999\n   68888\n')
f.close()

# Declaration du format
format = FortranFormat('2I4')

# Lecture 
import numpy as N
f  = open('misc.io.fortran.dat')
data = N.zeros((2, 2), 'i')
for i, line in enumerate(f):
    data[i] = FortranLine(line[:-1], format)[:]
print data
#  -> [[   5 9999]
#  ->  [   6 8888]]
