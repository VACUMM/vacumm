#!/usr/bin/env python
# -*- coding: utf8 -*-
"""De fortran à python avec f2py"""


# Importation et aide
import mymodule
help(mymodule)
print mymodule.myfunc.__doc__


# Utilisation
import numpy as N
a = N.arange(3.)
print mymodule.myfunc(a)
