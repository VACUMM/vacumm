# -*- coding: utf8 -*-
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
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

# NAR data / FTP CERSAT .... until 2009


class Old_nar(Sst) :
    """ NAR SST """
    def __init__(self):
        import os
        Sst.__init__(self)
        DIR_NAR=os.path.join(self.WORKDIR,'OLD_SST_NAR')      # repertoire de stockage des donnees SST NAR
        if os.path.isdir(DIR_NAR)==False:
            os.mkdir(DIR_NAR)
        self.WORKDIR=os.path.join(self.WORKDIR,'OLD_SST_NAR')

        #print "Workdir in nar.py"
        #print self.WORKDIR
        #ZONE_NAR='mocc' # Mediterranee Occidentale a relier a la congif MENOR
        self.ZONE_NAR='mocc'
        #self.ZONE_NAR='gasc' # Gascogne a relier a la config MANGA
        #ZONE_NAR='mnord' # Mer du Nord a relier a la config MANGA

        #print self.ZONE_NAR

        os.chdir(self.WORKDIR)

    def rappatrie(self,cfg=None):
        """ Rappatrie par ftp les donnees NAR SST """
        import os,subprocess,cdtime
        #----------------------------------------------------
        print ''
        print '---------- RECUPERATION FICHIERS NAR ----------'
        print ''
        #----------------------------------------------------


        #-------------------------------------------------------------
        #------- recuperation des donnees : generalites (site ftp, etc)
        if cfg is None:
	        config = ConfigParser.RawConfigParser()
	        config.read(os.path.join(self.SCRIPT_DIR,'config.cfg'))
	        URL_CERSAT = config.get('Old Nar SST', 'url_cersat')
	        #URL_CERSAT='ftp.ifremer.fr/pub/ifremer/cersat'
	        DATA_DIR = config.get('Old Nar SST', 'data_dir')
	        URL_NAR_DATA=os.path.join(URL_CERSAT,DATA_DIR)
	        GRID_DIR = config.get('Old Nar SST', 'grid_dir')
	        URL_NAR_GRID=os.path.join(URL_CERSAT,GRID_DIR)
	        znar=self.ZONE_NAR
	        gridfile="grid_%(znar)s.nc"%vars()
	        COMPRESS_NAR = config.get('Old Nar SST', 'compress')
	        EXT_NAR="00.nc.%(COMPRESS_NAR)s"%vars()
	        HEURE_NAR = config.get('Old Nar SST', 'nar_hours')
	        HEURE_NAR = HEURE_NAR.split(',')
	        NAR_DAYNIGHT = config.get('Old Nar SST', 'nar_daynight')    # NIGHT, DAY ou ALL
        else:
	        URL_CERSAT = cfg['Old Nar SST']['url_cersat']
	        #URL_CERSAT='ftp.ifremer.fr/pub/ifremer/cersat'
	        DATA_DIR = cfg['Old Nar SST']['data_dir']
	        URL_NAR_DATA=os.path.join(URL_CERSAT,DATA_DIR)
	        GRID_DIR = cfg['Old Nar SST']['grid_dir']
	        URL_NAR_GRID=os.path.join(URL_CERSAT,GRID_DIR)
	        znar=self.ZONE_NAR
	        gridfile="grid_%(znar)s.nc"%vars()
	        COMPRESS_NAR = cfg['Old Nar SST']['compress']
	        EXT_NAR="00.nc.%(COMPRESS_NAR)s"%vars()
	        HEURE_NAR = cfg['Old Nar SST']['nar_hours']
	        HEURE_NAR = HEURE_NAR.split(',')
	        NAR_DAYNIGHT = cfg['Old Nar SST']['nar_daynight']    # NIGHT, DAY ou ALL



        #-- recuperation de la grille NAR
        #----------------------------------------------------
        if os.path.isfile(gridfile)==False:
          url_dir="ftp://%(URL_NAR_GRID)s"%vars()
          print "Downloading NAR %(znar)s grid FILE"%vars()
          file_in="%(gridfile)s.gz"%vars()
          remote_file="%(url_dir)s/%(gridfile)s.gz"%vars()
          subprocess.call(["wget", "-O", file_in, remote_file])
          subprocess.call(["gunzip","-f", file_in])
        else:
            dwork=self.WORKDIR
            print "Grille NAR deja presente dans le repertoire %(dwork)s."%vars()


        #-- recuperation des donnees par FTP anonymous
        #----------------------------------------------------
        ctest=self.ctdeb

        # prevoir un test pour les cas ou le fichier de donnees n'existe pas !!!
        while ctest <= self.ctfin:
          if (NAR_DAYNIGHT == 'NIGHT') or (NAR_DAYNIGHT == 'ALL'):
            print '--> donnees de nuit'
            #H1=HEURE_NAR[0]
            H2=HEURE_NAR[3]
            #Old_nar.recup_ftp_sst_nar(self,URL_NAR_DATA,self.ZONE_NAR,ctest.year,ctest.month,ctest.day,H1, EXT_NAR, '9')
            Old_nar.recup_ftp_sst_nar(self,URL_NAR_DATA,self.ZONE_NAR,ctest.year,ctest.month,ctest.day,H2, EXT_NAR, '9')

          if (NAR_DAYNIGHT == 'DAY') or (NAR_DAYNIGHT == 'ALL'):
            print '--> donnees de jour'
            H1=HEURE_NAR[1]
            H2=HEURE_NAR[2]
            Old_nar.recup_ftp_sst_nar(self,URL_NAR_DATA,self.ZONE_NAR,ctest.year,ctest.month,ctest.day,H1, EXT_NAR, '9')
            Old_nar.recup_ftp_sst_nar(self,URL_NAR_DATA,self.ZONE_NAR,ctest.year,ctest.month,ctest.day,H2, EXT_NAR, '9')

          # On incremente le temps "test" de 1 jour
          ctest=ctest.add(1,cdtime.Days)

        #-- Fin de recuperation des donnees
        #----------------------------------------------------

    def recup_ftp_sst_nar(self,URL_NAR_DATA,ZONE,YYYY,MM,DD,HH,ext,ntry):
        """ Recuperation par FTP d un fichier de SST NAR horaire"""
        #**********************************************************************************
        # MOD : recup_ftp_sst_nar
        # OBJ : Recuperation par FTP d un fichier de SST NAR horaire
        # CRE : P. Craneguy (Actimar) / Projet PRECOC (ANR)
        #       G. Charria (Ifremer)
        # VER : 1.0 (25/01/2007)
        #       2.0 (11/2009)
        #**********************************************************************************
        import os,subprocess

        if MM < 10:
          MM="0%(MM)s"%vars()

        if DD < 10:
          DD="0%(DD)s"%vars()

        #-- URL site et Nom fichier a recuperer
        url_dir_nar="ftp://%(URL_NAR_DATA)s"%vars()
        url_dir="%(url_dir_nar)s/%(ZONE)s/%(YYYY)s"%vars()
        url_file_def="%(YYYY)s%(MM)s%(DD)s%(HH)s%(ext)s"%vars()
        dwork=self.WORKDIR
        outfile="%(dwork)s/%(url_file_def)s"%vars()

        #-- Extraction du fichier
        if  os.path.isfile(outfile)==False:
          print "Downloading NAR %(ZONE)s FILE %(url_file_def)s"%vars()
          list_file="%(url_dir)s/%(url_file_def)s"%vars()
          #subprocess.call(["wget", "-t", ntry, "-O", outfile, list_file])
          subprocess.call(["wget", "-t", ntry, "-O", outfile, list_file])
          s=os.path.getsize(outfile) # If file is empty !!
          if s==0:
            print "Empty file: %(outfile)s !"%vars()
            os.remove(outfile)
          else:
            subprocess.call(["gunzip","-f", outfile])

        else:
          print "Pas de downloading : %(url_dir)s/%(url_file_def)s existe deja"%vars()
        #-- Fin de ftp
        #----------------------------------------------------

    def read(self):
        """ Lecture des fichiers NetCDF de NAR SST """
        import cdms2,sys,os, glob
        import numpy,MV2
        import cdtime
        from vacumm.misc.axes import create_lon
        from vacumm.misc.grid import create_grid,  set_grid
        from vacumm.misc.atime import create_time
        from vacumm.misc.phys.units import kel2degc

        # -- Dans le cas d'un NAR SST on lit la grille lon lat dans le fichier grid
        # -- Lecture de la grille
        znar=self.ZONE_NAR
        gridfile="grid_%(znar)s.nc"%vars()
        f=cdms2.open(gridfile)
        la = f('latitude')
        lo = f('longitude')
        f.close()

        # -- Creation des axes longitude et latitude
        lat_axis = create_lon(la,id='latitude')
        lon_axis = create_lon(lo,id='longitude')
        # -- Creation de la grille
        grid = create_grid(lon_axis, lat_axis)

        # -- Creation d'un objet cdms nomme self.data et d'un tableau cumt pour les dates extraites des noms de fichiers
        self.data = () #Initialise un tuple
        # =============== ATTENTION ====================
        # Initialiser self.data pour ne pas dupliquer en memoire !!!!!!!!!
        # ============================================

        #cumt = [] # Initialise un array

        # -- Boucle sur les fichiers presents dans le WORKDIR

        #url_file_def="%(YYYY)s%(MM)s%(DD)s%(HH)s%(ext)s"%vars()

        #self.ctdeb

        # -- Ancienne methode
        #files = glob.glob(os.path.join(self.WORKDIR, '2*.nc'))
        #files.sort()
        # --

        # -- Methode amelioree
        # Cree une liste
        files = []
        # Cree la liste des fichiers correspondants a la periode consideree
        ctest = self.ctdeb
        while ctest <= self.ctfin:
            flnme_only = '%(#)04d%(##)02d%(###)02d*.nc'%{'#':ctest.year, '##':ctest.month, '###':ctest.day}
            files.extend(glob.glob(os.path.join(self.WORKDIR, flnme_only)))
            ctest=ctest.add(1,cdtime.Days)
        # --

        for filename in files:
            # -- Lecture du fichier filename
            f = cdms2.open(filename)
            temp = f('sea_surface_temperature')
            # =============== ATTENTION ==================================
            # Verifier que temp utilise tout le temps le meme espace memoire ... voir du cote de cdms
            # ==========================================================
            f.close()

            # -- Extraction de la date et heure du nom du fichier
            #ty = numpy.int(filename[-15:-11])
            #tm = numpy.int(filename[-11:-9])
            #tj = numpy.int(filename[-9:-7])
            #th = numpy.int(filename[-7:-5])
            #tt = cdtime.comptime(ty,tm,tj,th)
            #cumt.append(tt)

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



        # -- Attribution de la grille a l'objet self.data
        set_grid(self.data, grid, axes=True)

        # -- Informations sur le dataset
        self.data.name = "NAR_SST"
        self.data.units = "degree_Celsius"
        self.data.standard_name = "satellite_sea_surface_temperature"
        self.data.long_name = "Satellite Sea Surface Temperature - NAR"

        # -- Change unit
        self.data = kel2degc(self.data)

        #-- Fin de lecture des donnees
        #----------------------------------------------------
