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



class Seviri(Sst) :
    """ SEVIRI SST """
    def __init__(self,cfg=None):
        import os
        print 'using custom Seviri class !!'
        Sst.__init__(self,cfg)
        DIR_SEVIRI=os.path.join(self.WORKDIR,'SST_SEVIRI')      # repertoire de stockage des donnees SST SEVIRI
        if os.path.isdir(DIR_SEVIRI)==False:
            os.mkdir(DIR_SEVIRI)
        self.WORKDIR=os.path.join(self.WORKDIR,'SST_SEVIRI')
        self.OBSDIR=self.WORKDIR

        os.chdir(self.WORKDIR)
            
        

    def rappatrie(self,cfg=None):
        """ Rappatrie par ftp les donnees SEVIRI SST """
        import os,subprocess,cdtime
        #----------------------------------------------------
        print ''
        print '---------- RECUPERATION FICHIERS SEVIRI ----------'
        print ''
        #----------------------------------------------------
	
	#-------------------------------------------------------------
	#------- recuperation des donnees : generalites (site ftp, etc)
        if cfg is None:
  	    config = ConfigParser.RawConfigParser()
	    config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))
	
	try:
          if cfg is None: 
	      timerange = config.get('Seviri SST', 'timerange')
          else:
              timerange = cfg['Seviri SST']['timerange']
	  # timerange = 'midnight' pour donnees a minuit seulement 
	  # timerange = 'all' pour donnees a minuit seulement
	except ConfigParser.NoOptionError:
	  print 'No Time Range'
	  timerange = 'all' # Par defaut, lecture de toutes les heures
	
	if self.ctdeb >= cdtime.comptime(2012,01,01,0,0,0):
	  print timerange
	  print 'Data moved to /home5/taveeg/cache/project/osi-saf/data/sst/l3c/seviri/osi-saf/'
	  EXT_SEVIRI=".nc"
	  
	  if cfg is None:
	      URL_SEVIRI_DATA = config.get('Seviri SST', 'recent_dir')
	      URL_CERSAT = config.get('Seviri SST', 'url_cersat_rt')        
	      DATA_DIR = config.get('Seviri SST', 'data_dir_rt')
	      #URL_SEVIRI_DATA=os.path.join(URL_CERSAT,DATA_DIR)
	      
	      usr = config.get('Seviri SST', 'user') 
	      pwd = config.get('Seviri SST', 'pwd') 
          else:              
              URL_SEVIRI_DATA = cfg['Seviri SST']['recent_dir']
              URL_CERSAT = cfg['Seviri SST']['url_cersat_rt']
              DATA_DIR = cfg['Seviri SST']['data_dir_rt']
              #URL_SEVIRI_DATA=os.path.join(URL_CERSAT,DATA_DIR)

              usr = cfg['Seviri SST']['user']
              pwd = cfg['Seviri SST']['pwd']
	  
	  #-- recuperation des donnees par NFS
	  #----------------------------------------------------
	  ctest=self.ctdeb
	 
	  while ctest <= self.ctfin:   
	    # Pour le FTP fichiers telecharges corrompus ... pb BZ2
	    #Seviri.recup_ftp_sst_seviri_recent(self,URL_CERSAT,DATA_DIR,ctest.year,ctest.month,ctest.day, EXT_SEVIRI, usr,pwd, '9')
	    Seviri.recup_nfs_sst_seviri_recent(self,URL_SEVIRI_DATA,ctest.year,ctest.month,ctest.day,timerange)	    
	    
	    # On incremente le temps "test" de 1 jour
	    ctest=ctest.add(1,cdtime.Days)
	  
	else:
	  print 'Works for data before 2012/03/19. [In this program, method not used for data after 2012/01/01]'

          if cfg is None:
	      COMPRESS_SEVIRI = config.get('Seviri SST', 'compress')
	      URL_CERSAT = config.get('Seviri SST', 'url_cersat')        
	      DATA_DIR = config.get('Seviri SST', 'data_dir')
          else:
              COMPRESS_SEVIRI = cfg['Seviri SST']['compress']
              URL_CERSAT = cfg['Seviri SST']['url_cersat']
              DATA_DIR = cfg['Seviri SST']['data_dir']

	  URL_SEVIRI_DATA=os.path.join(URL_CERSAT,DATA_DIR)
	  EXT_SEVIRI=".nc.%(COMPRESS_SEVIRI)s"%vars()	
          
	  #NAR_DAYNIGHT = config.get('Nar SST', 'nar_daynight')    # NIGHT, DAY ou ALL
	  
	  #data_dir = products/gridded/experimental-cms/netcdf/msg
	  #url_cersat = ftp.ifremer.fr/pub/ifremer/cersat
	  #compress = bz2
	  
      
      

	  #-- recuperation des donnees par FTP anonymous
	  #----------------------------------------------------
	  ctest=self.ctdeb

	  # prevoir un test pour les cas ou le fichier de donnees n'existe pas !!! 
	  while ctest <= self.ctfin:  
	    Seviri.recup_ftp_sst_seviri(self,URL_SEVIRI_DATA,ctest.year,ctest.month,ctest.day, EXT_SEVIRI, '9',timerange,cfg)
	    
	    # On incremente le temps "test" de 1 jour
	    ctest=ctest.add(1,cdtime.Days)
        
       
        #-- Fin de recuperation des donnees
        #----------------------------------------------------
        
    def recup_nfs_sst_seviri_recent(self,URL_SEVIRI_DATA,YYYY,MM,DD,timerange):
	""" Recuperation par NFS d un fichier de SST SEVIRI horaire recent (depuis 01/01/2012)""" 
	import os, subprocess,shutil
	import cdtime,glob
	from vacumm.misc.atime import strftime
	
	# -- Complement du nom de repertoire avec l'annee et le jour relatif (001-> 365)
	DATA_DIR = os.path.join(URL_SEVIRI_DATA,str(YYYY))
	#print DATA_DIR
	a = cdtime.comptime(YYYY,MM,DD)
	a2 = a.torel('days since '+strftime('%Y-%m-%d',cdtime.comptime(YYYY-1,12,31)))
	
	DATA_DIR = os.path.join(DATA_DIR,'%(#)03d' % {'#':a2.value})
	
	#print DATA_DIR
	if timerange == 'midnight':
	  filename = '%(#)04d%(##)02d%(####)02d00*' % {'#':YYYY, '##':MM, '####':DD}    
	else:	  
	  filename = '%(#)04d%(##)02d%(####)02d*' % {'#':YYYY, '##':MM, '####':DD}
	copy_mode = 'nfs'
	list_file = glob.glob(os.path.join(DATA_DIR,filename))
	# si la liste est vide, on essaie en se connectant a service7
	if not list_file:
	  find_cmd = 'ssh service7 "find %(DATA_DIR)s -name \'%(filename)s\' "' %vars()
	  list_file = subprocess.check_output(find_cmd, shell=True).strip().split()
	  if list_file:
	    copy_mode = 'scp'
	
	for file_to_read in list_file:	
	  if os.path.isfile(os.path.basename(file_to_read)) == False :
	    if copy_mode == 'nfs':
	      shutil.copyfile(file_to_read, os.path.basename(file_to_read))
	    if copy_mode == 'scp':
	      copy_cmd = 'scp caparmor-sftp.ifremer.fr:%(file_to_read)s .' %vars()
	      subprocess.check_call( copy_cmd, shell=True)
        
        #-- Fin de ftp
        #----------------------------------------------------
        
    def recup_ftp_sst_seviri_recent(self,URL_CERSAT,DATA_DIR,YYYY,MM,DD,ext,usr,pwd,ntry,timerange):
	""" Recuperation par FTP d un fichier de SST SEVIRI horaire recent (depuis 01/01/2012)""" 
	import os, subprocess
	from ftplib import FTP
	import cdtime
	from vacumm.misc.atime import strftime
	
	
	
	# --------------- LE FTP FONCTIONNE MAIS LES FICHIERS SEMBLENT CORROMPUS ... PROBLEME DANS LE BZ2
	# -- Complement du nom de repertoire avec l'annee et le jour relatif (001-> 365)
	DATA_DIR = os.path.join(DATA_DIR,str(YYYY))
	print DATA_DIR
	a = cdtime.comptime(YYYY,MM,DD)
	a2 = a.torel('days since '+strftime('%Y-%m-%d',cdtime.comptime(YYYY-1,12,31)))
	
	DATA_DIR = os.path.join(DATA_DIR,'%(#)03d' % {'#':a2.value})
	
	print DATA_DIR
	
	# connection au CDOCO
        ftp = FTP(host=URL_CERSAT, user=usr, passwd=pwd)
        
        # positionnement sur le repertoire "best_estimate"
        ftp.cwd(DATA_DIR)
        
        # utile pour ftp.retrbinary
        def handleDownload(block):
            file.write(block)
	
	if timerange == 'midnight':
	  filename = '%(#)04d%(##)02d%(####)02d00*' % {'#':YYYY, '##':MM, '####':DD}    
	else:	  
	  filename = '%(#)04d%(##)02d%(####)02d*' % {'#':YYYY, '##':MM, '####':DD}
        #print filename
        
        list_file = ftp.nlst(filename)
        #print list_file
        
        for file_to_read in list_file:
            # recuperation fichier si non present dans work/MODEL/MARS
            if os.path.isfile(file_to_read) == False:
                # rajouter test sur date du fichier ... telecharge que si plus recent ...
                
                file = open(file_to_read, 'wb')
                
                ftp.retrbinary('RETR ' + file_to_read , handleDownload)
                
                subprocess.call(["bunzip2","-f", file_to_read])
	
	    else:
	      print "Pas de downloading : %(file_to_read)s existe deja"%vars()
        
        ftp.quit()
        # ----------------------------------------------------------------------------------------
        
        #-- Fin de ftp
        #----------------------------------------------------
	

    def recup_ftp_sst_seviri(self,URL_SEVIRI_DATA,YYYY,MM,DD,ext,ntry,timerange,cfg):
        """ Recuperation par FTP d un fichier de SST SEVIRI horaire""" 
        #**********************************************************************************
        # MOD : recup_ftp_sst_seviri
        # OBJ : Recuperation par FTP d un fichier de SST SEVIRI horaire
        # CRE : G. Charria (Ifremer)
        # VER : 1.0 (11/10/2010)
        #**********************************************************************************
        import os,subprocess
        
        if MM < 10:
          MM="0%(MM)s"%vars()
       
        if DD < 10:
          DD="0%(DD)s"%vars()
        
        #-- URL site et Nom fichier a recuperer
        url_dir_seviri="ftp://%(URL_SEVIRI_DATA)s"%vars()
        url_dir="%(url_dir_seviri)s/%(YYYY)s/%(MM)s"%vars()
        url_file_def="sst1h_msg_%(YYYY)s%(MM)s%(DD)s%(ext)s"%vars()
        if cfg is None:
            config = ConfigParser.RawConfigParser()
            config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))
            COMPRESS_SEVIRI = config.get('Seviri SST', 'compress')
        else:
            COMPRESS_SEVIRI = cfg['Seviri SST']['compress']
        extu = ext[:-(len(COMPRESS_SEVIRI)+1)]
        url_file_defu="sst1h_msg_%(YYYY)s%(MM)s%(DD)s%(extu)s"%vars()
        dwork=self.WORKDIR
        outfile="%(dwork)s/%(url_file_def)s"%vars()
        outfile_uncompressed="%(dwork)s/%(url_file_defu)s"%vars()
        
        # S. Petton
        # On regarde si on souhaite restreindre l'emprise geographique
        if 'latmin' in cfg['Observations']:
            latmin = cfg['Observations']['latmin']
        else:
            latmin=None
        if 'latmax' in cfg['Observations']:
            latmax = cfg['Observations']['latmax']
        else:
            latmax=None
        if 'lonmin' in cfg['Observations']:
            lonmin = cfg['Observations']['lonmin']
        else:
            lonmin=None
        if 'lonmax' in cfg['Observations']:
            lonmax = cfg['Observations']['lonmax']
        else:
            lonmax=None
        # S. Petton
 
        #-- Extraction du fichier
        if  os.path.isfile(outfile_uncompressed)==False:
          print "Downloading SEVIRI FILE %(url_file_def)s"%vars()
          list_file="%(url_dir)s/%(url_file_def)s"%vars()
          #subprocess.call(["wget", "-t", ntry, "-O", outfile, list_file])
          subprocess.call(["wget", "-t", ntry, "-O", outfile, list_file])
          s=os.path.getsize(outfile) # If file is empty !!
          if s==0:
            print "Empty file: %(outfile)s !"%vars()
            os.remove(outfile)
          else:
            subprocess.call(["bunzip2","-f", outfile])
            # S. Petton
            if lonmin != None and lonmax != None and latmin != None and latmax != None :
               outncfile=outfile[:len(outfile)-4]
               nco_cmd = "ncks -h -a -O %s -o tmp.nc" %(outncfile)  
               nco_cmd +=  " -d lon,%s,%s -d lat,%s,%s "%(lonmin,lonmax,latmin,latmax)
               subprocess.call(nco_cmd,shell=True,stdin=None,stdout=None)
               mv_cmd = "mv tmp.nc %s" %(outncfile)
               subprocess.call( mv_cmd,shell=True,stdin=None,stdout=None)
            # S. Petton
    
        else:
          print "Pas de downloading : %(url_dir)s/%(url_file_def)s existe deja"%vars()
        #-- Fin de ftp
        #----------------------------------------------------

    def read(self,cfg=None):
        """ Lecture des fichiers NetCDF de NAR SST """
        import cdms2,sys,os, glob
        import numpy,MV2
        import cdtime
        from vacumm.misc.axes import create_lon
        from vacumm.misc.grid import create_grid,  set_grid
        from vacumm.misc.atime import create_time
        from vacumm.misc.phys.units import kel2degc
        
        if self.ctdeb >= cdtime.comptime(2012,01,01,0,0,0):
	  # -- Creation d'un objet cdms nomme self.data et d'un tableau cumt pour les dates extraites des noms de fichiers
	  self.data = () #Initialise un tuple
	  # =============== ATTENTION ====================
	  # Initialiser self.data pour ne pas dupliquer en memoire !!!!!!!!!
	  # ============================================
	  
	  # -- Methode amelioree
	  # Cree une liste
	  files = []
	  # Cree la liste des fichiers correspondants a la periode consideree
	  ctest = self.ctdeb
	  while ctest <= self.ctfin:
	      for iH in numpy.arange(24): 
		flnme_only = '%04d%02d%02d%02d*.nc'%(ctest.year, ctest.month, ctest.day, iH)
		files.extend(glob.glob(os.path.join(self.OBSDIR, flnme_only)))
	      ctest=ctest.add(1,cdtime.Days)
	    # --

          if cfg is None:
	      config = ConfigParser.RawConfigParser()
	      config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))
	      lomin = float(config.get('Domain', 'lomin') )    
	      lomax = float(config.get('Domain', 'lomax')      )
	      lamin = float(config.get('Domain', 'lamin')     )
	      lamax = float(config.get('Domain', 'lamax')     )    
          else:
              lomin = cfg['Domain']['lomin']
              lomax = cfg['Domain']['lomax']
              lamin = cfg['Domain']['lamin']
              lamax = cfg['Domain']['lamax']

	  if cfg is None:
	    try: 
	      timerange = config.get('Seviri SST', 'timerange')
	      # timerange = 'midnight' pour donnees a minuit seulement 
	      # timerange = 'all' pour donnees a minuit seulement
	    except ConfigParser.NoOptionError:
	      #print 'No Time Range'
	      timerange = 'all' # Par defaut, lecture de toutes les heures
	  else:
	    timerange = cfg['Seviri SST']['timerange']
	  
	  if files == []:
	      print 'No data file to read ...'
	  else:        
	      for ifile, filename in enumerate(files):
		  # -- Lecture du fichier filename                        
		  f = cdms2.open(filename)
		  temp = f('sea_surface_temperature', lon=(lomin,lomax), lat=(lamin,lamax))		      
		  f.close()
		  
		  
		  # -- Transfert de temp dans l'objet cdat self.data (concatenation)		  
		  self.data += temp, 
		  
	      
	      # -- Creation de l'axe temporel 
	      #taxis = create_time(cumt)
	      
	      # -- MV2.concatenate pour concatener les informations dans self.data (entre autre construit l'axe temporel)
	      self.data = MV2.concatenate(self.data)
					      
	      # -- Informations sur le dataset
	      #self.data.name = "SEVIRI_SST"
	      self.data.units = "degree_Celsius"
	      self.data.standard_name = "satellite_sea_surface_temperature"
	      self.data.long_name = "Satellite Sea Surface Temperature - SEVIRI"
	      
	      # -- Change unit
	      self.data = kel2degc(self.data)
	    
        else:
	  # -- Dans le cas d'un NAR SST on lit la grille lon lat dans le fichier grid
	  
	  # -- Creation d'un objet cdms nomme self.data et d'un tableau cumt pour les dates extraites des noms de fichiers
	  self.data = () #Initialise un tuple
	  # =============== ATTENTION ====================
	  # Initialiser self.data pour ne pas dupliquer en memoire !!!!!!!!!
	  # ============================================
	  
	  # -- Methode amelioree
	  # Cree une liste
	  files = []
	  # Cree la liste des fichiers correspondants a la periode consideree
	  ctest = self.ctdeb
	  while ctest <= self.ctfin:
	      flnme_only = '*%04d%02d%02d*.nc'%(ctest.year, ctest.month, ctest.day)
	      files.extend(glob.glob(os.path.join(self.OBSDIR, flnme_only)))
	      ctest=ctest.add(1,cdtime.Days)
	  # --
          
          if cfg is None:
              config = ConfigParser.RawConfigParser()
              config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))
              lomin = float(config.get('Domain', 'lomin') )
              lomax = float(config.get('Domain', 'lomax')      )
              lamin = float(config.get('Domain', 'lamin')     )
              lamax = float(config.get('Domain', 'lamax')     )
          else:
              lomin = cfg['Domain']['lomin']
              lomax = cfg['Domain']['lomax']
              lamin = cfg['Domain']['lamin']
              lamax = cfg['Domain']['lamax']

          if cfg is None:
	    try: 
	      timerange = config.get('Seviri SST', 'timerange')
	      # timerange = 'midnight' pour donnees a minuit seulement 
	      # timerange = 'all' pour donnees a minuit seulement
	    except ConfigParser.NoOptionError:
	      #print 'No Time Range'
	      timerange = 'all' # Par defaut, lecture de toutes les heures
	  else:
	    timerange = cfg['Seviri SST']['timerange']
	
	  
	  if files == []:
	      print 'No data file to read ...'
	  else:        
	      for ifile, filename in enumerate(files):
		  # -- Lecture du fichier filename                        
		  f = cdms2.open(filename)
		  
		  if timerange=='midnight':
		    temp = f('sst', lon=(lomin,lomax), lat=(lamin,lamax), time=slice(0,1))   
		  else:
		    temp = f('sst', lon=(lomin,lomax), lat=(lamin,lamax))
		  
		  # =============== ATTENTION ==================================
		  # Verifier que temp utilise tout le temps le meme espace memoire ... voir du cote de cdms
		  # ==========================================================
		  f.close()
		  
		  
		  # -- Transfert de temp dans l'objet cdat self.data (concatenation)
		  # =============== ATTENTION ====================
		  # Faire une variable intermediaire qui sera au prealable allouee en memoire pour eviter
		  # trop de copies en memoire !!!!!!!!!!!!!!!!
		  # ============================================
		  self.data += temp, 
		  
	      
	      # -- Creation de l'axe temporel 
	      #taxis = create_time(cumt)
	      
	      # -- MV2.concatenate pour concatener les informations dans self.data (entre autre construit l'axe temporel)
	      self.data = MV2.concatenate(self.data)
					      
	      # -- Informations sur le dataset
	      #self.data.name = "SEVIRI_SST"
	      self.data.units = "degree_Celsius"
	      self.data.standard_name = "satellite_sea_surface_temperature"
	      self.data.long_name = "Satellite Sea Surface Temperature - SEVIRI"
	      
	      # -- Change unit
	      self.data = kel2degc(self.data)
        
        #-- Fin de lecture des donnees
        #----------------------------------------------------
        
    def read_assim(self,cfg=None):
        """ Lecture des fichiers NetCDF de SST SEVIRI formatte pour l'assimilation """
        import cdms2,sys,os, glob
        import numpy,MV2
        import cdtime
                
        # -- Creation d'un objet cdms nomme self.data                 
                
        if cfg is None:
	  config = ConfigParser.RawConfigParser()
	  config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))
	  lomin = float(config.get('Domain', 'lomin') )    
	  lomax = float(config.get('Domain', 'lomax')      )
	  lamin = float(config.get('Domain', 'lamin')     )
	  lamax = float(config.get('Domain', 'lamax')     )   
	  obsdir = config.get('Observations',  'obsdir')
	  fln = config.get('Observations',  'filename')
	else:
	  lomin = cfg['Domain']['lomin']
	  lomax = cfg['Domain']['lomax']
	  lamin = cfg['Domain']['lamin']
	  lamax = cfg['Domain']['lamax']
	  obsdir = cfg['Observations']['obsdir']
	  fln = cfg['Observations']['filename']
	  
        
        filename = os.path.join(obsdir, fln) 
        
        if os.path.isfile(filename) == False :                    
            print 'No data file to read ...'
        else:        
            # -- Lecture du fichier filename                        
            f = cdms2.open(filename)
            self.data = f('temp', lon=(lomin,lomax), lat=(lamin,lamax),  time=(self.ctdeb, self.ctfin))     
            f.close()                        
                                            
        # -- Informations sur le dataset
        #self.data.name = "SEVIRI_SST"
        self.data.units = "degree_Celsius"
        self.data.standard_name = "satellite_sea_surface_temperature"
        self.data.long_name = "Satellite Sea Surface Temperature - SEVIRI"            
    def read_long(self,cfg=None):
        """ Lecture des fichiers NetCDF de SST SEVIRI formatte pour runs longs """
        import cdms2,sys,os, glob
        import numpy,MV2
        import cdtime
        from vacumm.misc.io import list_forecast_files
        from vacumm.misc.phys.units import kel2degc
        
        lomin = cfg['Domain']['lomin']
        lomax = cfg['Domain']['lomax']
        lamin = cfg['Domain']['lamin']
        lamax = cfg['Domain']['lamax']
        obsdir = cfg['Observations']['obsdir']
        fln = cfg['Observations']['filename']
        
        # -- Lecture du fichier filename
        self.data=[]
        flist = list_forecast_files(os.path.join(obsdir,fln),(str(cfg['Time Period']['andeb'])+'-'+str(cfg['Time Period']['mdeb'])+'-'+str(cfg['Time Period']['jdeb']),str(cfg['Time Period']['anfin'])+'-'+str(cfg['Time Period']['mfin'])+'-'+str(cfg['Time Period']['jfin']),'co'))  
        print flist
        for fil in flist:                      
            f = cdms2.open(fil)
            temp = f('sst', lon=(lomin,lomax), lat=(lamin,lamax),  time=(self.ctdeb, self.ctfin))     
            f.close()
            
            if len(flist) > 1:
                self.data += temp,   
            else:
                self.data = temp          
            
        if len(flist) > 1:
            # -- MV2.concatenate pour concatener les informations dans self.data (entre autre construit l'axe temporel)
            self.data = MV2.concatenate(self.data)
        
        
        # -- Change unit
        self.data = kel2degc(self.data)
            
        # -- Informations sur le dataset
        #self.data.name = "SEVIRI_SST"
        self.data.units = "degree_Celsius"
        self.data.standard_name = "satellite_sea_surface_temperature"
        self.data.long_name = "Satellite Sea Surface Temperature - SEVIRI"            
        
        #-- Fin de lecture des donnees
        #----------------------------------------------------
        
        
        
