# -*- coding: utf8 -*-
"""


"""
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

#from ftplib import FTP

# utile pour ftp.retrbinary
def handleDownload(block):
    file.write(block)
    #print ".",

def get_ftp_f1(ctdeb, ctfin,andeb, SCRIPT_DIR, cfg=None):
    import os,cdtime, subprocess, ConfigParser
    print 'ftp cdoco en cours'
    # connection au F1/F2 au CDOCO

    if cfg is None:
      config = ConfigParser.RawConfigParser()
      config.read(os.path.join(SCRIPT_DIR,'config.cfg'))

      if config.get('Model Description', 'name') == 'mars_manga':
        dirf1 = config.get('MARS F1', 'url_f1')
        usr = config.get('MARS F1', 'user')
        pwd = config.get('MARS F1', 'pwd')
        fic_prefix = config.get('MARS F1', 'fic_prefix')       #'PREVIMER_F1-MARS3D-MANGA4000_'
        time_res = float(config.get('MARS F1', 'time_res') )   # pas de temps des sorties en heure

      if config.get('Model Description', 'name') == 'mars_menor':
        dirf1 = config.get('MARS F2', 'url_f2')
        usr = config.get('MARS F2', 'user')
        pwd = config.get('MARS F2', 'pwd')
        fic_prefix = config.get('MARS F2', 'fic_prefix')       #'PREVIMER_F2-MARS3D-MENOR_'
        time_res = float(config.get('MARS F2', 'time_res') )   # pas de temps des sorties en heure

    else:

      if cfg['Model Description']['name'] == 'mars_manga':
        dirf1 = cfg['MARS F1']['url_f1']
        usr = cfg['MARS F1']['user']
        pwd = cfg['MARS F1']['pwd']
        fic_prefix = cfg['MARS F1']['fic_prefix']       #'PREVIMER_F1-MARS3D-MANGA4000_'
        time_res = cfg['MARS F1']['time_res']   # pas de temps des sorties en heure

      if cfg['Model Description']['name'] == 'mars_menor':
        dirf1 = cfg['MARS F2']['url_f2']
        usr = cfg['MARS F2']['user']
        pwd = cfg['MARS F2']['pwd']
        fic_prefix = cfg['MARS F2']['fic_prefix']       #'PREVIMER_F2-MARS3D-MENOR_'
        time_res = cfg['MARS F2']['time_res']   # pas de temps des sorties en heure




    #ftp = FTP(host='eftp.ifremer.fr',user='c1e975',passwd='Ef0XmnZ4')
    # positionnement sur le repertoire "best_estimate"
    #ftp.cwd('./f1_4000/best_estimate/')
    # positionnement sur le repertoire de l'annee
    #ftp.cwd(str(andeb))

    base_filename = os.path.join(dirf1,str(andeb))

    ctest = ctdeb
    current_year = ctdeb.year

    while ctest <= ctfin:

        # si on change d'annee alors on change de repertoire
        if ctest.year != current_year:
            #ftp.cwd(str(ctest.year))
            current_year = ctest.year
            base_filename = os.path.join(dirf1,str(ctest.year))

        # ajoute un zero si hdeb > 10 pour obtenir la chaine '03' au lieu de '3' par exemple :
        #if ctest.hour < 10:
        #    hour='0%s'%(ctest.hour)
        #else:
        #    hour=str(ctest.hour)


        # a faire aussi pour mois et jour !!!

        #construction du nom de fichier
        fic = fic_prefix+'%4d%02d%02dT%02d00Z.nc' %(ctest.year,ctest.month,ctest.day,ctest.hour)
        filename = ''
        filename=os.path.join(base_filename,fic)
        #print filename

        # recuperation fichier si non present dans work/MODEL/MARS
        if os.path.isfile(fic)==False:

            full_filename = 'ftp://%(usr)s:%(pwd)s@%(filename)s' %vars()


            print full_filename
            subprocess.call(["wget", full_filename])

        else:
            print "No ftp copy needed : %(fic)s already exists."%vars()

        # Open the file for writing in binary mode
        #file = open(filename, 'wb')

        #ftp.retrbinary('RETR ' +filename,handleDownload)
        # executer un test pour voir si le fichier existe dans la base de donnee
        # ....
        # ....


        # On incremente le temps "test" de 1 ou  3 heures (time_res)
        ctest=ctest.add(time_res,cdtime.Hours)

    #ftp.quit()
