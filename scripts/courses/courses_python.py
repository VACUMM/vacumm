#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Python de base"""

# Importations
import numpy
a = numpy.arange(5)
from numpy import arange
b = arange(5)


# Inspection
help(numpy.arange)
print dir(numpy)
print type('3')
print callable(numpy.arange)
print isinstance('a', str)


# Indentations
a = True
if a is True:
    print 'True !!!'        # -> ESSAYER SANS INDENTATION


# Quelques types

# - listes
a = [1, 3] 
a.append(5)                 # -> ESSAYER EXTEND()
a[0] = 10
print a
del a[1]

# - tuples
b = (5, 6)
print b[1]                  # -> ESSAYER DE CHANGER B[0], DEL
c = list(a)
print b, c

# - dictionnaires
c = {'a':1,  'b':2}
d = dict(c=4)
c.update(d, e=5)
print c.keys()
print c


# Boucles

# - for
a = range(6)
for i in a:
    print i
for key, val in c.items():
    print key, val
# - while
i = 0
while i < 5:
    i += 1
print i


# Interception d'erreur
try:
    print AAAAA
except:
    print 'Erreur'
    
    
# IMPORTANT : REFERENCE AUX VARIABLES NON PRIMAIRES

# - lien
a = [3, 5]
b = a
print a is b
b[0] = 100
print a

# - copy
ac = a.copy()
from copy import copy, deepcopy
ac2 = copy(a)
adc = deepcopy(a)


# Fonctions

# - inline
myfunc = lambda arg: arg+1
print myfunc(2)

# - normale
def myfunc2(arg, kwarg=4):
    """Documentation"""
    print arg, kwarg
# -> TESTER HELP, KWARG
print myfunc2.__doc__

# - arguments speciaux
def myfunc3(*args, **kwargs):
    print 'args', args
    print 'kwargs', kwargs
def myfunc4(a, b=1, **kwargs):
    print a, b, kwargs
# -> TESTER


# Classes

# - declaration
class MaClasse(object):
    """Doc"""
    def __init__(self, a):
        self.a = a
    def foo(self, b):
        """Doc foo"""
        return self.a+b
        
# - initialisation
mc = MaClasse(5)
print type(mc)
print isinstance(mc, MaClasse)

# - utilisation
print mc.foo(10), mc.foo(100)


# Destruction
del mc
