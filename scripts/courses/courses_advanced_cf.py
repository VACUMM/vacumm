#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Utilitaires et conventions de formatage

"""

from vcmq import MV2, cdms2, N
from vacumm.data.cf import VAR_SPECS, format_var, match_var, format_axis, AXIS_SPECS, GENERIC_VAR_NAMES


# Création d'une variable 1D
sst = MV2.arange(5.)                        # -> VERIFIER LES INFOS DE CET VARIABLE
sst.getAxis(0).designateLongitude()         # -> VERIFIER LES INFOS DE CET AXE


# Formatage
format_var(sst, 'sst')                      # -> VERIFIER LES NOUVELLES INFOS DE VAR+AXE
# -> FORMATER L'AXE EN LATITUDE AU POINT U


# Verification
print match_var(sst,'sst')
