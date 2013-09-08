#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Passer par un logger"""

from vacumm.misc.log import Logger # logger avancé
#from vacumm.misc.io import Logger # logger plus simple
from vcmq import dict_filter


# Initialisation
logger = Logger('TEST', level='info', 
    format='[%(asctime)s %(name)s %(levelname)s] %(message)s')
    
    
# Test
logger.info('mon info')
logger.warning('mon warning')
logger.debug('mon debug')                   # -> RETESTER AVEC UN SET_LEVEL POUR DEBUG
logger.error('mon erreur')

        


