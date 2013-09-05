#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Gestion de configurations avancées"""

from vcmq import ConfigManager
from ConfigParser import SafeConfigParser


# Initialisation avec fichier de specifications
cfgm = ConfigManager('courses_advanced_cfgm.ini')


# Valeurs par défaut
print cfgm.defaults()                           # -> CHANGER LE .INI


# Chargement de sa configuration personnelle
cfg = cfgm.load('courses_advanced_cfgm.cfg')    # -> CHANGER LE .CFG
print cfg
print type(cfg)


# Modification
cfg['plot']['linewidth'] = 2.
print cfg


# -> AJOUTEZ UNE SOUS-SECTION AUX .INI ET .CFG


# Via ConfigParser
cfg2 = SafeConfigParser()
cfg2.read('courses_advanced_cfgm.cfg')
print cfg2 
print cfg2.get('plot', 'linewidth')
print cfg2.getfloat('plot', 'linewidth')

