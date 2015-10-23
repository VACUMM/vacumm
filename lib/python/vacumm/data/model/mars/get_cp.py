
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

import sys, os, shutil
import cdtime

def get_cp(ctdeb, ctfin, dir_f1, cfg=None):

    if cfg is None:
        print 'Not implemented for a use without cfg/ini files.'
    sys.exit()

    ctest = ctdeb

    while ctest <= ctfin:

        # si on change d'annee alors on change de repertoire
        #if ctest.year != current_year:
        #   ftp.cwd(str(ctest.year))
        #  current_year = ctest.year

        # ajoute un zero si hdeb > 10 pour obtenir la chaine '03' au lieu de '3' par exemple :
        if ctest.hour < 10:
            hour='0%s'%(ctest.hour)
        else:
            hour=str(ctest.hour)

         # a faire aussi pour mois et jour

        #construction du nom de fichier
        #filename='PREVIMER_F1-MARS3D-MANGA_'+str(ctest.year)+str(ctest.month)+str(ctest.day)+'T'+hour+'00Z.nc'
        #A changer pour que le chemin d'accès ne soit plus en dur en reutilisable pour Mars Menor. A lire direct dans config.gfg.
        if cfg['Model Description']['name'] == 'mars_manga':
            filename=cfg['MARS F1']['fic_prefix']+'%4d%02d%02dT%02d00Z.nc' %(ctest.year,ctest.month,ctest.day,ctest.hour)
        if cfg['Model Description']['name'] == 'mars_menor':
            filename=cfg['MARS F2']['fic_prefix']+'%4d%02d%02dT%02d00Z.nc' %(ctest.year,ctest.month,ctest.day,ctest.hour)
        print filename

        if os.path.isfile(filename) == False :
            shutil.copyfile(dir_f1+filename, filename)
            print filename, 'doesn t exists'


    # On incremente le temps "test" de 1 heure
        if cfg['Model Description']['name'] == 'mars_manga':
                inc=cfg['MARS F1']['time_res']
        if cfg['Model Description']['name'] == 'mars_menor':
                inc=cfg['MARS F2']['time_res']

        ctest=ctest.add(inc,cdtime.Hours)

def get_cp_f1(ctdeb, ctfin, dir_f1):
    ctest = ctdeb

    while ctest <= ctfin:

        # si on change d'annee alors on change de repertoire
        #if ctest.year != current_year:
        #   ftp.cwd(str(ctest.year))
        #  current_year = ctest.year

        # ajoute un zero si hdeb > 10 pour obtenir la chaine '03' au lieu de '3' par exemple :
        if ctest.hour < 10:
            hour='0%s'%(ctest.hour)
        else:
            hour=str(ctest.hour)

         # a faire aussi pour mois et jour

        #construction du nom de fichier
        #filename='PREVIMER_F1-MARS3D-MANGA_'+str(ctest.year)+str(ctest.month)+str(ctest.day)+'T'+hour+'00Z.nc'
        #A changer pour que le chemin d'accès ne soit plus en dur en reutilisable pour Mars Menor. A lire direct dans config.gfg.
        filename='PREVIMER_F1-MARS3D-MANGA4000_%4d%02d%02dT%02d00Z.nc' %(ctest.year,ctest.month,ctest.day,ctest.hour)
        print filename

        if os.path.isfile(filename) == False :
            shutil.copyfile(dir_f1+filename, filename)
            print filename, 'doesn t exists'


    # On incremente le temps "test" de 1 heure
        ctest=ctest.add(1,cdtime.Hours)


def get_cp_f2_v9(ctdeb, ctfin, dir_f2_v9):
    ctest = ctdeb

    while ctest <= ctfin:

        # si on change d'annee alors on change de repertoire
        #if ctest.year != current_year:
        #   ftp.cwd(str(ctest.year))
        #  current_year = ctest.year

        # ajoute un zero si hdeb > 10 pour obtenir la chaine '03' au lieu de '3' par exemple :
        if ctest.hour < 10:
            hour='0%s'%(ctest.hour)
        else:
            hour=str(ctest.hour)

         # a faire aussi pour mois et jour

        #construction du nom de fichier
        #filename='PREVIMER_F1-MARS3D-MANGA_'+str(ctest.year)+str(ctest.month)+str(ctest.day)+'T'+hour+'00Z.nc'
        #A changer pour que le chemin d'accès ne soit plus en dur en reutilisable pour Mars Menor. A lire direct dans config.gfg.
        filename='PREVIMER_F2-MARS3D-MENOR1200_%4d%02d%02dT%02d00Z.nc' %(ctest.year,ctest.month,ctest.day,ctest.hour)
        print filename

        if os.path.isfile(filename) == False :
            shutil.copyfile(dir_f2_v9+filename, filename)


    # On incremente le temps "test" de 3 heure
        ctest=ctest.add(3,cdtime.Hours)
