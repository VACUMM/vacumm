# Classe satellite: rappatriement (dans le WORKDIR) et lecture des observations satellites
#**********************************************************************************
# CRE : P. Craneguy (Actimar) - Projet PRECOC (ANR)
#        G. Charria (Ifremer)
#        S. Theetten (Ifremer)
# VER : 1.0 (25/01/2007)
#        2.0 (11/2009) - Python
#       3.0 (03/2010) - passage en classes
#**********************************************************************************
# --------------------------------------------
## Erreurs lorsque on active ces deux lignes.
## => UnboundLocalError: local variable 'parent' referenced before assignment
# --------------------------------------------
#from vacumm.misc.bases import Object

#class Satellite(Object):
# --------------------------------------------
class Satellite():
    """ Donnees Satellite """
    def __init__(self):
        print ''
