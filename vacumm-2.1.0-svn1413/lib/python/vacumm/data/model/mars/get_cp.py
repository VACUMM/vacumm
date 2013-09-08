
# -*- coding: utf8 -*-
"""
    

"""
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
