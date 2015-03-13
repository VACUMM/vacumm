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
from vacumm.data.satellite.sst import Sst
import ConfigParser

#import psutil

# NAR data / Recent data ... not on ftp

class Nar(Sst) :
    """ NAR SST """
    def __init__(self):
        import os
        Sst.__init__(self)
        DIR_NAR=os.path.join(self.WORKDIR,'SST_NAR')      # repertoire de stockage des donnees SST NAR  
        if os.path.isdir(DIR_NAR)==False:
            os.mkdir(DIR_NAR)
        self.WORKDIR=os.path.join(self.WORKDIR,'SST_NAR')
            
        os.chdir(self.WORKDIR)

    def rappatrie(self,cfg=None):
        """ Rappatrie par cp des donnees NAR SST AVHRR METOP_A"""
        import os,subprocess,cdtime
        #----------------------------------------------------
        print '---------- RECUPERATION FICHIERS NAR SST ----------'
        #----------------------------------------------------
        
        if cfg is None:
	  config = ConfigParser.RawConfigParser()
	  config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))        
   
    
        

        #-- recuperation des donnees 
        #--------------------------------------
        ctest=self.ctdeb

        # prevoir un test pour les cas ou le fichier de donnees n'existe pas !!! 
        while ctest <= self.ctfin:     
        
	    if cfg is None:
        
	      if config.get('Observations', 'download_detail') == 'nfs':
		  Nar.recup_cp_sst(self,ctest)
	      if config.get('Observations', 'download_detail') == 'ftp':  
		  Nar.recup_ftp_sst(self,ctest)
		  
	    else:
        
	      if cfg['Observations']['download_detail'] == 'nfs':
		  Nar.recup_cp_sst(self,ctest)
	      if cfg['Observations']['download_detail'] == 'ftp':  
		  Nar.recup_ftp_sst(self,ctest)
	      
            # On incremente le temps "test" de 1 jour
            ctest=ctest.add(1,cdtime.Days)
       
        #-- Fin de recuperation des donnees
        #----------------------------------------------------

    def recup_cp_sst(self,ctps,cfg=None):
        """ Recuperation par COPY d un fichier de SST NAR de Nuit""" 
        import os,  shutil,  subprocess
              
        #-------------------------------------------------------------
        #------- recuperation des donnees : generalites (repertoire, etc)
        if cfg is None:
	  config = ConfigParser.RawConfigParser()
	  config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))        
	  file_directory = config.get('Nar SST', 'dir_nar')
	else:       
	  file_directory = cfg['Nar SST']['dir_nar']
	  
        
        #-- Repertoire et Nom fichier a recuperer        
        dwork=self.WORKDIR
        julian_day = ctps.torel('days since %4d'%(ctps.year))
        julian_day.value = julian_day.value+1
        
        if cfg is None:
	  satellites = config.get('Nar SST', 'satellites')
	  satellites = satellites.split(',')
	  satellites_detail = config.get('Nar SST', 'satellites_detail')
	  satellites_detail = satellites_detail.split(',')
	  satellites_dir = config.get('Nar SST', 'satellites_dir')
	  satellites_dir = satellites_dir.split(',')
	  hr_satellites = config.get('Nar SST', 'hr_satellites')
	  hr_satellites = hr_satellites.split(',')
	else:
	  satellites = cfg['Nar SST']['satellites']
	  satellites = satellites.split(',')
	  satellites_detail = cfg['Nar SST']['satellites_detail']
	  satellites_detail = satellites_detail.split(',')
	  satellites_dir = cfg['Nar SST']['satellites_dir']
	  satellites_dir = satellites_dir.split(',')
	  hr_satellites = cfg['Nar SST']['hr_satellites']
	  hr_satellites = hr_satellites.split(',')
	  
        for isat, s_name in enumerate(satellites):
            url_dir="%s/%s/%4d/%03d/"%(file_directory, satellites_dir[isat], ctps.year, julian_day.value)
            file_def="%4d%02d%02d-EUR-L3P_GHRSST-SSTsubskin-NAR_AVHRR_%s-eumetsat_sstnar_%s_%4d%02d%02d_%s-v01.7-fv01.0.nc.bz2" %(ctps.year, 
                                                                                                                                 ctps.month, ctps.day,s_name, satellites_detail[isat],   ctps.year, ctps.month, ctps.day,  hr_satellites[isat])                
            #print file_def
            # Exemples de noms de fichiers:
            #20100121-EUR-L3P_GHRSST-SSTsubskin-NAR_AVHRR_NOAA_19-eumetsat_sstnar_noaa19_20100121_130000-v01.7-fv01.0.nc.bz2            
            #20100121-EUR-L3P_GHRSST-SSTsubskin-NAR_AVHRR_METOP_A-eumetsat_sstnar_metop02_20100121_200000-v01.7-fv01.0.nc.bz2  
            outfile="%(file_def)s"%vars()            
            
            #-- Extraction du fichier
            # outfile[:-4] => filename without bz2 extension
            if  os.path.isfile(outfile[:-4])==False:    
                print url_dir+file_def        
                if os.path.exists(url_dir+file_def):
                    print "Copying NAR FILE %(file_def)s"%vars()            
                    shutil.copyfile(url_dir+file_def, file_def)                
                    subprocess.call(["bunzip2","-f", outfile])                
                else:
                    print "The file does not exist ! => %(url_dir)s%(file_def)s"%vars()
            else:
                print "No copy needed : %(url_dir)s%(file_def)s already exists."%vars()
                
        #-- Fin de copie
        #----------------------------------------------------
        
    def recup_ftp_sst(self,ctps,cfg=None):
        """ Recuperation par FTP d un fichier de SST NAR de Nuit""" 
        import os,  shutil,  subprocess
              
        #-------------------------------------------------------------
        #------- recuperation des donnees : generalites (repertoire, etc)
        if cfg is None:
	  config = ConfigParser.RawConfigParser()
	  config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))        
	  URL_CERSAT = config.get('Nar SST', 'url_cersat')
	  DATA_DIR = config.get('Nar SST', 'data_dir')
	  file_directory = URL_CERSAT+DATA_DIR     # file_directory = config.get('Nar SST', 'dir_nar')
	  usr = config.get('Nar SST', 'user')
	  pwd = config.get('Nar SST', 'pwd')
	  full_file_directory = 'ftp://%(usr)s:%(pwd)s@%(file_directory)s' %vars()
	  
	  #-- Repertoire et Nom fichier a recuperer        
	  julian_day = ctps.torel('days since %4d'%(ctps.year))
	  julian_day.value = julian_day.value+1
	  
	  satellites = config.get('Nar SST', 'satellites')
	  satellites = satellites.split(',')
	  satellites_detail = config.get('Nar SST', 'satellites_detail')
	  satellites_detail = satellites_detail.split(',')
	  satellites_dir = config.get('Nar SST', 'satellites_dir')
	  satellites_dir = satellites_dir.split(',')
	  hr_satellites = config.get('Nar SST', 'hr_satellites')
	  hr_satellites = hr_satellites.split(',')
	else:      
	  URL_CERSAT = cfg['Nar SST']['url_cersat']
	  DATA_DIR = cfg['Nar SST']['data_dir']
	  file_directory = URL_CERSAT+DATA_DIR     # file_directory = config.get('Nar SST', 'dir_nar')
	  usr = cfg['Nar SST']['user']
	  pwd = cfg['Nar SST']['pwd']
	  full_file_directory = 'ftp://%(usr)s:%(pwd)s@%(file_directory)s' %vars()
	  
	  #-- Repertoire et Nom fichier a recuperer        
	  julian_day = ctps.torel('days since %4d'%(ctps.year))
	  julian_day.value = julian_day.value+1
	  
	  satellites = cfg['Nar SST']['satellites']
	  satellites = satellites.split(',')
	  satellites_detail = cfg['Nar SST']['satellites_detail']
	  satellites_detail = satellites_detail.split(',')
	  satellites_dir = cfg['Nar SST']['satellites_dir']
	  satellites_dir = satellites_dir.split(',')
	  hr_satellites = cfg['Nar SST']['hr_satellites']
	  hr_satellites = hr_satellites.split(',')
	  
        for isat, s_name in enumerate(satellites):
            url_dir="%s/%s/%4d/%03d/"%(full_file_directory, satellites_dir[isat], ctps.year, julian_day.value)
            file_def="%4d%02d%02d-EUR-L3P_GHRSST-SSTsubskin-NAR_AVHRR_%s-eumetsat_sstnar_%s_%4d%02d%02d_%s-v01.7-fv01.0.nc.bz2" %(ctps.year, 
                                                                                                                                 ctps.month, ctps.day,s_name, satellites_detail[isat],   ctps.year, ctps.month, ctps.day,  hr_satellites[isat])                
            #print file_def
            # Exemples de noms de fichiers:
            #20100121-EUR-L3P_GHRSST-SSTsubskin-NAR_AVHRR_NOAA_19-eumetsat_sstnar_noaa19_20100121_130000-v01.7-fv01.0.nc.bz2            
            #20100121-EUR-L3P_GHRSST-SSTsubskin-NAR_AVHRR_METOP_A-eumetsat_sstnar_metop02_20100121_200000-v01.7-fv01.0.nc.bz2  
            outfile="%(file_def)s"%vars()            
            
            #-- Extraction du fichier
            # outfile[:-4] => filename without bz2 extension
            if  os.path.isfile(outfile[:-4])==False:    
                print url_dir+file_def        
                #if os.path.exists(url_dir+file_def):
                print "Downloading NAR FILE %(file_def)s"%vars() 
                subprocess.call(["wget", url_dir+file_def])
                subprocess.call(["bunzip2","-f", outfile])                
                #else:
                #    print "The file does not exist ! => %(url_dir)s%(file_def)s"%vars()
            else:
                print "No ftp copy needed : %(file_def)s already exists."%vars()
                
        #-- Fin de copie
        #----------------------------------------------------
        
    def read(self, verbose=False, cfg=None):
        """ Lecture des fichiers NetCDF de NAR SST """
        import cdms2,sys,os, glob
        import numpy,MV2
        import cdtime
        from vacumm.misc.axes import create_lon, set_order
        from vacumm.misc.grid import create_grid,  set_grid
        from vacumm.misc.atime import create_time
        from vacumm.misc.phys.units import kel2degc
        import gc
      
        # Get the configuration file information 
        #cfg = self.get_config()
        #print cfg
              
        # -- Creation d'un objet cdms nomme self.data et d'un tableau cumt pour les dates extraites des noms de fichiers
        self.data = () #Initialise un tuple
        # =============== ATTENTION ====================
        # Initialiser self.data pour ne pas dupliquer en memoire !!!!!!!!!
        # ============================================

        # -- Methode amelioree
        # Cree une liste
        files = []
        # Cree la liste des fichiers correspondants a la periode consideree
        if cfg is None:
	  config = ConfigParser.RawConfigParser()
	  config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))
	  hr_satellites = config.get('Nar SST', 'hr_satellites')
	else:
	  hr_satellites = cfg['Nar SST']['hr_satellites']
	  
	hr_satellites = hr_satellites.split(',')
	#hr_satellites = cfg['hr_satellites']
	  

        #print hr_satellites
        
        ctest = self.ctdeb
        while ctest <= self.ctfin:
            for isat, s_name in enumerate(hr_satellites):
                
                flnme_only = '%04d%02d%02d*%s*.nc'%(ctest.year, ctest.month, ctest.day, s_name)
                files.extend(glob.glob(os.path.join(flnme_only)))
                
            ctest=ctest.add(1,cdtime.Days)
        # --
      
	if cfg is None:
	  lomin = float(config.get('Domain', 'lomin') )
	  lomax = float(config.get('Domain', 'lomax')      )
	  lamin = float(config.get('Domain', 'lamin')     )
	  lamax = float(config.get('Domain', 'lamax')     )
	else:
	  lomin = cfg['Domain']['lomin']
	  lomax = cfg['Domain']['lomax']
	  lamin = cfg['Domain']['lamin']
	  lamax = cfg['Domain']['lamax']
	  
        
        #print files

        if files == []:
            print 'No data file to read ...'
        else:
            
            # ---- Lecture et creation de la grille ----
            #
            # -- Lecture du fichier filename
            f = cdms2.open(files[0])                        
            lo = f.getVariable('lon')
            la = f.getVariable('lat')
            # -- Creation des axes longitude et latitude
#            lat_axis = create_lon(la,id='latitude')
#            lon_axis = create_lon(lo,id='longitude')
            # -- Creation de la grille
            grid = create_grid(lo, la)
            
            del lo,  la
            
            for ifile, filename in enumerate(files):
      
                
                # -- Lecture du fichier filename
                f = cdms2.open(filename)
                
                temp2 = f('sea_surface_temperature')
                set_order(temp2, 'tyx') # pour que averager.py fontionne
                
                # modif J.Gatti : utilisation de l'index de qualite (0:unprocessed 1:cloudy 2:bad 3:suspect 4:acceptable 5:excellent)
                temp2.set_fill_value(temp2._FillValue)
                conf=f('proximity_confidence')
                MV2.putmask(temp2.data,conf.data<3,temp2._FillValue)   #  ne change que data   
                MV2.putmask(temp2.mask,conf.data<3,True)                    # ne change que mask
                #autre methode
                #--------------------
                #temp2=MV2.masked_where(conf.data<3,temp2)          # ne change que mask
                #oldmask=temp2.mask
                #temp2[:]=temp2.filled()                                                  # change data mais met mask a false partout
                #temp2.mask=oldmask                                                     # remet le bon mask sans lien
                del conf
                # fin modif J.Gatti : utilisation de l'index de qualite    
                
                # -- Attribution de la grille a l'objet self.data
                set_grid(temp2, grid, axes=True)     
                
                temp = temp2(lon=(lomin, lomax), lat=(lamin, lamax))
                
                if verbose:
                    # == TEST OCCUPATION MEMOIRE ===
                    print ctest,  'Avant'
                    #print psutil.Process(os.getpid()).get_memory_percent()
                    #print psutil.Process(os.getpid()).get_memory_info()
                    #print 'CPU percent: ', psutil.cpu_percent(interval=0.1)
                    #print 'Used phymem: ', psutil.used_phymem()
                    #print 'Used virtmem: ', psutil.used_virtmem()
                
                del temp2
                
                if verbose:
                    # == TEST OCCUPATION MEMOIRE ===
                    print ctest,  'Apres del'
                    #print psutil.Process(os.getpid()).get_memory_percent()
                    #print psutil.Process(os.getpid()).get_memory_info()
                    #print 'CPU percent: ', psutil.cpu_percent(interval=0.1)
                    #print 'Used phymem: ', psutil.used_phymem()
                    #print 'Used virtmem: ', psutil.used_virtmem()
 
                # =============== ATTENTION ==================================
                # Verifier que temp utilise tout le temps le meme espace memoire ... voir du cote de cdms
                # ==========================================================
                f.close()
                
                
                


                # -- Transfert de temp dans l'objet cdat self.data (concatenation)
                # =============== ATTENTION ====================
                # Faire une variable intermediaire qui sera au prealable allouee en memoire pour eviter
                # trop de copies en memoire !!!!!!!!!!!!!!!!
                # ============================================
                #self.data += temp,
                if ifile == 0:
                    self.data = temp
                else:                
                    self.data = MV2.concatenate((self.data, temp))
                
                if verbose:
                    # == TEST OCCUPATION MEMOIRE ===
                    print ctest,  'Avant gccollect'
                    #print psutil.Process(os.getpid()).get_memory_percent()
                    #print psutil.Process(os.getpid()).get_memory_info()
                    #print 'CPU percent: ', psutil.cpu_percent(interval=0.1)
                    #print 'Used phymem: ', psutil.used_phymem()
                    #print 'Used virtmem: ', psutil.used_virtmem()
                    print gc.collect()
                gc.collect()
                
                if verbose:
                    # == TEST OCCUPATION MEMOIRE ===
                    print ctest,  'Apres gccollect'
                    #print psutil.Process(os.getpid()).get_memory_percent()
                    #print psutil.Process(os.getpid()).get_memory_info()
                    #print 'CPU percent: ', psutil.cpu_percent(interval=0.1)
                    #print 'Used phymem: ', psutil.used_phymem()
                    #print 'Used virtmem: ', psutil.used_virtmem()

            # -- Creation de l'axe temporel
            #taxis = create_time(cumt)
            
            # -- MV2.concatenate pour concatener les informations dans self.data (entre autre construit l'axe temporel)
            #self.data = MV2.concatenate(self.data)
 
            
 
 
        # -- Informations sur le dataset
        #self.data.name = "NAR_SST"
        self.data.units = "degree_Celsius"
        self.data.standard_name = "satellite_sea_surface_temperature"
        self.data.long_name = "Satellite Sea Surface Temperature - NAR"
        
        # -- Change unit
        self.data = kel2degc(self.data)
        
        #-- Fin de lecture des donnees
        #----------------------------------------------------

    def clean_data(self):
        """ Delete data files to save disk space """
        import os, glob, cdtime

        # -- Methode amelioree
        # Cree une liste
        files = []
        # Cree la liste des fichiers correspondants a la periode consideree
        ctest = self.ctdeb
        while ctest <= self.ctfin:
            flnme_only = '%04d%02d%02d*.nc'%(ctest.year, ctest.month, ctest.day)
            files.extend(glob.glob(os.path.join(flnme_only)))
            ctest=ctest.add(1,cdtime.Days)
        # --

        if files == []:
            print 'No data file to delete ...'
        else:
            for filename in files:
                os.remove(filename)
        
        
        
        
    
