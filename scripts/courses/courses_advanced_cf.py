#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Utilitaires et conventions de formatage (:mod:`vacumm.data.cf`) :

- Définitions des prioriétés standards de variables usuelles comme une
  extension aux convention CF.
- Outils de formatage.


.. note:: Les spécifications seront probablement converties au format xml à l'avenir.

"""

from vcmq import MV2, cdms2, N
from vacumm.data.cf import var_specs, format_var, match_var, format_axis, axis_specs, generic_var_names


# Création d'une variable 1D
sst = MV2.arange(5.)                        # -> VERIFIER LES INFOS DE CET VARIABLE
sst.getAxis(0).designateLongitude()         # -> VERIFIER LES INFOS DE CET AXE


# Formatage
format_var(sst, 'sst')                      # -> VERIFIER LES NOUVELLES INFOS DE VAR+AXE
# -> FORMATER L'AXE EN LATITUDE AU POINT U


# Verification
print match_var(sst,'sst')
