#!/usr/bin/env python
# -*- coding: utf8 -*-
#
# Copyright or © or Copr. Actimar/IFREMER (2013-2015)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
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
