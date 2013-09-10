#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Python de base"""

# Importations
import numpy
a = numpy.arange(5)
from numpy import arange
b = arange(5)


# Inspection
# -> help(numpy.arange)
print dir(numpy)
print type('3')
print callable(numpy.arange)
print isinstance('a', str)


# Indentations
a = True
if a is True:
    print 'True !!!'        # -> ESSAYEZ SANS INDENTATION


# Quelques types

# - listes
a = [1, 3] 
a.append(5)                 # -> ESSAYEZ EXTEND()
a[0] = 10
print a
del a[1]

# - tuples
b = (5, 6)
print b[1]                  # -> ESSAYEZ DE CHANGER B[0], DEL
c = list(a)
print b, c

# - dictionnaires
c = {'a':1,  'b':2}
d = dict(c=4)
c.update(d, e=5)
print c.keys()
print c

# - chaines de caracteres
a = 'a'
b = "b"
c = """
e
bc
"""
print "a=%s b=%i"%(3, 4)
a=3; b=4
print "a=%(a)s b=%(b)i"%dict(a=a, b=b)  # -> ESSAYEZ AVEC VARS()
print 'jean dupont'.title()
print '123'.rjust(100)                  # -> ESSAYEZ .CENTER()
print 'a={} b={}'.format(a, b)
print 'a={0} b={1:f}'.format(a, b)
a = {'arg1':'ok', 'arg2':6}
print 'arg1 de a = {a[arg1]}, arg2 = {a[arg2]}'.format(a=a) # on ne passe que a

# Slicing
a = range(10)
print a[0], a[-1]
print a[-2:]
print a[::2]                # -> COMMENCER A 3
print a[slice(None, None, 2)]


# Boucles

# - for
a = range(6)
for i in a:
    print i
    if i==2:
         continue
    elif i==4:
         break
d = {'a':1,  'b':2}
for key, val in d.items():
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

# - copy (rarement utilisé)
from copy import copy, deepcopy
ac = copy(a)
ec2 = list(a)
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
myfunc4(1, **dict(b=3, c=5))
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
