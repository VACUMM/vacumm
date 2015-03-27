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
from vacumm.data.satellite.satellite import Satellite
import ConfigParser

class Sst(Satellite) :
    """ SST satellite """
    def __init__(self,cfg=None):
        import os, cdtime
        Satellite.__init__(self)

	SCRIPT_DIR=os.getcwd()
        self.SCRIPT_DIR = SCRIPT_DIR

        #-- retrocompatibilite
        if cfg is None:
            config = ConfigParser.RawConfigParser()
            config.read(os.path.join(SCRIPT_DIR,'config.cfg'))	
            andeb = config.getint('Time Period', 'andeb')
            anfin = config.getint('Time Period', 'anfin')		
            mdeb = config.getint('Time Period', 'mdeb')
            mfin = config.getint('Time Period', 'mfin')
            jdeb = config.getint('Time Period', 'jdeb')
            jfin = config.getint('Time Period', 'jfin')
            hdeb = config.getint('Time Period', 'hdeb')
            hfin = config.getint('Time Period', 'hfin')
            
            self.WORKDIR = config.get('Env', 'workdir')
        else:
            
            andeb = cfg['Time Period']['andeb']
            anfin = cfg['Time Period']['anfin']
            mdeb = cfg['Time Period']['mdeb']
            mfin = cfg['Time Period']['mfin']
            jdeb = cfg['Time Period']['jdeb']
            jfin = cfg['Time Period']['jfin']
            hdeb = cfg['Time Period']['hdeb']
            hfin = cfg['Time Period']['hfin']

            self.WORKDIR = cfg['Env']['workdir']

        #print "%(andeb)s-%(mdeb)s-%(jdeb)s %(hdeb)s"%vars()
        #print "%(anfin)s-%(mfin)s-%(jfin)s %(hfin)s"%vars()

        # Conversion en "Component Time"
        self.ctdeb=cdtime.comptime(andeb,mdeb,jdeb,hdeb,0,0)
        self.ctfin=cdtime.comptime(anfin,mfin,jfin,hfin,0,0)


        DIR_SST=os.path.join(self.WORKDIR,'SST')    # repertoire de stockage des donnees     
        if os.path.isdir(DIR_SST)==False:
            os.mkdir(DIR_SST)
        self.WORKDIR=os.path.join(self.WORKDIR,'SST')
            
        self.name="Satellite Sea surface Temperature"
        self.shortname="SST"
        self.units=" "

        #print "Workdir in sst.py"
        #print self.WORKDIR
