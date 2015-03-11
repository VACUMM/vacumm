# -*- coding: utf8 -*-
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

# Classe recopesca: rappatriement (dans le WORKDIR) et lecture des profiles RECOPESCA
#**********************************************************************************
# CRE : G. Charria (Ifremer)
#        S. Theetten (Ifremer)
# VER : 1.0 (12/10/2010)
#**********************************************************************************
import os, sys, ConfigParser
from vacumm.data.in_situ.profile import Profile

class Recopesca():

    qc_ok = "0111111"

    """ Donnees RECOPESCA """
    def __init__(self):
        import cdtime

        SCRIPT_DIR = os.getcwd()
        #print SCRIPT_DIR
        #os.chdir('../')
        #BASENAME = os.getcwd()
        #print BASENAME
        #os.chdir(SCRIPT_DIR) # back to the script_dir
        #self.WORKDIR = os.path.join(BASENAME, 'work')
        # repertoire de travail
        #print self.WORKDIR
        self.SCRIPT_DIR = SCRIPT_DIR

        # Lecture de la periode de validation
        config = ConfigParser.RawConfigParser()
        config.read(os.path.join(SCRIPT_DIR, 'config.cfg'))
        andeb = config.getint('Time Period', 'andeb')
        anfin = config.getint('Time Period', 'anfin')
        mdeb = config.getint('Time Period', 'mdeb')
        mfin = config.getint('Time Period', 'mfin')
        jdeb = config.getint('Time Period', 'jdeb')
        jfin = config.getint('Time Period', 'jfin')
        hdeb = config.getint('Time Period', 'hdeb')
        hfin = config.getint('Time Period', 'hfin')

        self.WORKDIR = config.get('Env', 'workdir')

        # Conversion en "Component Time"
        self.ctdeb = cdtime.comptime(andeb, mdeb, jdeb, hdeb, 0, 0)
        self.ctfin = cdtime.comptime(anfin, mfin, jfin, hfin, 0, 0)

        DIR_REC = os.path.join(self.WORKDIR, 'RECOPESCA')    # repertoire de stockage des donnees
        if os.path.isdir(DIR_REC) == False:
            os.chdir(self.WORKDIR)
            os.mkdir('RECOPESCA')
        self.WORKDIR = os.path.join(self.WORKDIR, 'RECOPESCA')

        self.name = "RECOPESCA data"
        self.shortname = "Recopesca"

        #dictionnaire contient les profils (cle : date)
        self.map_profiles = {}

    def rappatrie(self):
        """ Rappatrie par ftp les donnees RECOPESCA """
        import subprocess, cdtime

        os.chdir(self.WORKDIR)  # on se place dans le repertoire de travail
        #----------------------------------------------------
        print '---------- RECUPERATION FICHIERS RECOPESCA ----------'
        print 65 * '-'
        #----------------------------------------------------


        #-------------------------------------------------------------
        #------- recuperation des donnees : generalites (site ftp, etc)
        config = ConfigParser.RawConfigParser()
        config.read(os.path.join(self.SCRIPT_DIR, 'config.cfg'))
        URL_CDOCO = config.get('RECOPESCA', 'url_cdoco')
        DATA_DIR = config.get('RECOPESCA', 'data_dir')
        URL_REC_DATA = os.path.join(URL_CDOCO, DATA_DIR)
        EXT = ".csv"
        usr = config.get('RECOPESCA', 'user')
        pwd = config.get('RECOPESCA', 'pwd')

        os.chdir(self.WORKDIR)

        #-- recuperation des donnees par FTP anonymous
        #----------------------------------------------------
        ctest = self.ctdeb

        # prevoir un test pour les cas ou le fichier de donnees n'existe pas !!!
        while ctest.month <= self.ctfin.month:

            Recopesca.recup_ftp_recopesca(self, URL_CDOCO, DATA_DIR, ctest.year, ctest.month, EXT, usr, pwd, '9')
            # On incremente le temps "test" de 1 mois (va au premier du mois suivant)
            ctest = ctest.add(1, cdtime.Months)

        os.chdir(self.WORKDIR)
        #-- Fin de recuperation des donnees
        #----------------------------------------------------

    def recup_ftp_recopesca(self, URL_DATA, DATA_DIR, YYYY, MM, ext, usr, pwd, ntry):
        """ Recuperation par FTP d un fichier RECOPESCA mensuel"""
        from ftplib import FTP

        if MM < 10:
          MM = "0%(MM)s" % vars()

        # connection au CDOCO
        ftp = FTP(host=URL_DATA, user=usr, passwd=pwd)

        # positionnement sur le repertoire "best_estimate"
        ftp.cwd(DATA_DIR)

        # utile pour ftp.retrbinary
        def handleDownload(block):
            file.write(block)

        filename = '*%(#)s%(##)s*%(###)s' % {'#':YYYY, '##':MM, '###':ext}
        #print filename

        list_file = ftp.nlst(filename)
        #print list_file

        for file_to_read in list_file:
            # recuperation fichier si non present dans work/MODEL/MARS
            if os.path.isfile(file_to_read) == False:
                # rajouter test sur date du fichier ... telecharge que si plus recent ...

                file = open(file_to_read, 'wb')

                ftp.retrbinary('RETR ' + file_to_read , handleDownload)
            # executer un test pour voir si le fichier existe dans la base de donnee
            # ....
            # ....

        ftp.quit()


        #-- Fin de ftp
        #----------------------------------------------------

    def read(self):
        """ Lecture des fichiers csv de RECOPESCA """
        import glob
        #from vacumm.misc.atime import strptime
        from vacumm.misc.atime import numtime, strtime
        from matplotlib.dates import num2date
        import cdtime

        # - convertisseur pour le temps
        def convtime(s):
            return numtime(s)

        #----------------------------------------------------
        ctest = self.ctdeb
        ext = ".csv"
        platform_name = list()
        while ctest.month <= self.ctfin.month:

            if ctest.month < 10:
                MM = "0%(#)s" % {'#':ctest.month}
            else:
                MM = "%(#)s" % {'#':ctest.month}

            # On incremente le temps "test" de 1 mois
            filename = '*%(#)s%(##)s*%(###)s' % {'#':ctest.year, '##':MM, '###':ext}
            print glob.glob(filename)

            #creaton map temporaire
            map_profiles_tmp = {}

            #parcours des fichiers
            for file in glob.glob(filename):

                #lecture fichier
                file_data = open(file).read()
                #suppression ligne 1
                first_line, file_data = file_data.split('\n', 1)

                #extraction des profil du fichier (ligne a ligne)
                for line in file_data.splitlines():

                    line = line.strip().split(',')

                    platform_name = line[0]
                    time = line[1]
                    lat = line[2]
                    lon = line[3]
                    depth = line[4]
                    temp = line[5]
                    sal = line[6]
                    qc = line[7]

                    #test qc (si != 0111111 non traite)
                    if qc != self.qc_ok:
                        continue

                    time = num2date(convtime(time))
                    key = strtime(time)

                    #si profil existe
                    if key in map_profiles_tmp:
                        profile = map_profiles_tmp.get(key)
                        profile.add_time_step(depth, temp, sal)

                    #sinon creation profil
                    else:
                        profile = Profile(time, platform_name, lon, lat)
                        profile.add_time_step(depth, temp, sal)
                        map_profiles_tmp[key] = profile

                #creation tables numpy
                for key in map_profiles_tmp:
                    profile=map_profiles_tmp[key]
                    profile.convert_tables()

                #copie des profil du fichier
                self.map_profiles.update(map_profiles_tmp)

                #reinitialisation
                map_profiles_tmp.clear()

            ctest = ctest.add(1, cdtime.Months) # Passe au premier du mois suivant

    def get_profiles_position(self):
        from matplotlib.dates import date2num
        from numpy import array
        time_list = []
        lon_list = []
        lat_list = []

        #creation list
        for key in self.map_profiles:
            profile=self.map_profiles[key]
            time, lon, lat = profile.get_position()
            time_list.append(date2num(time))
            lat_list.append(lat)
            lon_list.append(lon)

        #conversion list table numpy
        time_list = array(time_list, dtype='float')
        lon_list = array(lon_list, dtype='float')
        lat_list = array(lat_list, dtype='float')

        return [time_list, lon_list, lat_list]










