
# Module de rappatriement des donnees SST dans le WORKDIR
#**********************************************************************************
# MOD : r_sst_nar.py
# OBJ : Recuperation par FTP anonymous des fichiers de SST NAR sur une 
#        periode de 'ndayfile' jours depuis la date 'DATE_DEB' puis concatenation de
#       ces fichiers dans 'FILESORT_NAR'
#
# CRE : P. Craneguy (Actimar) - Projet PRECOC (ANR)
#        G. Charria (Ifremer)
#        S. Theetten (Ifremer)
# VER : 1.0 (25/01/2007)
#        2.0 (11/2009) - Python
#**********************************************************************************
def rappatrie(WORKDIR,DIR_SST_NAR,andeb,mdeb,jdeb,hdeb,anfin,mfin,jfin,hfin,ZONE_NAR):
    import os,subprocess
    from cdtime import s2c
    import cdtime
    #import cdtime
    os.chdir(WORKDIR)  # on se place dans le repertoire de travail
    #----------------------------------------------------
    print ''
    print '---------- RECUPERATION FICHIERS NAR ----------'
    print ''
    #----------------------------------------------------


    #-------------------------------------------------------------
    #------- recuperation des donnees : generalites (site ftp, etc)
    URL_CERSAT='ftp.ifremer.fr/pub/ifremer/cersat'
    URL_NAR_DATA=os.path.join(URL_CERSAT,'SAFOSI/Products/NARSST/netcdf')
    URL_NAR_GRID=os.path.join(URL_CERSAT,'SAFOSI/Products/NARSST/netcdf/grids')
    gridfile="grid_%(ZONE_NAR)s.nc"%vars()
    COMPRESS_NAR='gz'
    EXT_NAR="00.nc.%(COMPRESS_NAR)s"%vars()
    HEURE_NAR=['02','10','12','20']
    NAR_DAYNIGHT='NIGHT'	# NIGHT, DAY ou ALL



    os.chdir(DIR_SST_NAR)
    #-- recuperation de la grille NAR
    #----------------------------------------------------
    if os.path.isfile(gridfile)==False:
      url_dir="ftp://%(URL_NAR_GRID)s"%vars()
      print "Downloading NAR %(ZONE_NAR)s grid FILE"%vars()
      file_in="%(gridfile)s.gz"%vars()
      remote_file="%(url_dir)s/%(gridfile)s.gz"%vars()
      subprocess.call(["wget", "-O", file_in, remote_file])
      subprocess.call(["gunzip","-f", file_in])
    else:
      print "Grille NAR deja presente dans le repertoire %(DIR_SST_NAR)s."%vars()


    #-- recuperation des donnees par FTP anonymous
    #----------------------------------------------------

    # Conversion en "Component Time"
    ctdeb=s2c("%(andeb)s-%(mdeb)s-%(jdeb)s %(hdeb)s"%vars())
    ctfin=s2c("%(anfin)s-%(mfin)s-%(jfin)s %(hfin)s"%vars())

    ctest=ctdeb

    # prevoir un test pour les cas ou le fichier de donnees n'existe pas !!! 
    while ctest <= ctfin:
      print ctest
      if (NAR_DAYNIGHT == 'NIGHT') or (NAR_DAYNIGHT == 'ALL'):
        print '--> donnees de nuit'
        H1=HEURE_NAR[0]
        H2=HEURE_NAR[3]
        recup_ftp_sst_nar(URL_NAR_DATA,ZONE_NAR,ctest.year,ctest.month,ctest.day,H1, EXT_NAR, DIR_SST_NAR,'9')
        recup_ftp_sst_nar(URL_NAR_DATA,ZONE_NAR,ctest.year,ctest.month,ctest.day,H2, EXT_NAR, DIR_SST_NAR,'9')

      if (NAR_DAYNIGHT == 'DAY') or (NAR_DAYNIGHT == 'ALL'): 
        print '--> donnees de jour'
        H1=HEURE_NAR[1]
        H2=HEURE_NAR[2]      
        recup_ftp_sst_nar(URL_NAR_DATA,ZONE_NAR,ctest.year,ctest.month,ctest.day,H1, EXT_NAR, DIR_SST_NAR,'9')
        recup_ftp_sst_nar(URL_NAR_DATA,ZONE_NAR,ctest.year,ctest.month,ctest.day,H2, EXT_NAR, DIR_SST_NAR,'9')
      
      # On incremente le temps "test" de 1 jour
      ctest=ctest.add(1,cdtime.Days)
   
    os.chdir(WORKDIR)
    #-- Fin de recuperation des donnees
    #----------------------------------------------------

def recup_ftp_sst_nar(URL_NAR_DATA,ZONE,YYYY,MM,DD,HH,ext,DIR_SST_NAR,ntry):
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
    print DD

    if DD < 10:
      DD="0%(DD)s"%vars()

    print DD

    #-- URL site et Nom fichier a recuperer
    url_dir_nar="ftp://%(URL_NAR_DATA)s"%vars()
    url_dir="%(url_dir_nar)s/%(ZONE)s/%(YYYY)s"%vars()
    url_file_def="%(YYYY)s%(MM)s%(DD)s%(HH)s%(ext)s"%vars()
    outfile="%(DIR_SST_NAR)s/%(url_file_def)s"%vars()
    
    #-- Extraction du fichier
    if  os.path.isfile(outfile)==False:
      print "Downloading NAR %(ZONE)s FILE %(url_file_def)s"%vars()
      list_file="%(url_dir)s/%(url_file_def)s"%vars()
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



def read(WORKDIR,DIR_SST_NAR,andeb,mdeb,jdeb,hdeb,anfin,mfin,jfin,hfin,ZONE_NAR):
    import cdms2,sys,os,vcs
    import numpy,MV2
    import cdtime
 
    os.chdir(DIR_SST_NAR)
    # -- Lecture de la grille
    gridfile="grid_%(ZONE_NAR)s.nc"%vars()
    f=cdms2.open(gridfile)
    lat_obs = f('latitude')
    lon_obs = f('longitude')
    f.close()
    
    # -- Lecture de la SST
    list_file=os.listdir(DIR_SST_NAR)
    list_file.sort()
    icount=0

    # Maladroit mais pas d'autre solution pour le moment
    # pour estime la taille du vecteur temps
    #for xfile in list_file:
    #  if (andeb in xfile) or (anfin in xfile):
    #    icount=icount+1
    #time_obs=numpy.empty((icount))
    #y=numpy.empty((icount))
    #m=numpy.empty((icount))
    #j=numpy.empty((icount))
    #h=numpy.empty((icount))


    icount=0
    sst_sat=[]
    qc_sst_sat=[]
    time_obs=[]
    for xfile in list_file:
      if (andeb in xfile) or (anfin in xfile):
        print "Fichier lu: %(xfile)s"%vars()
        f=cdms2.open(xfile)
        sst_temp = f('sea_surface_temperature',squeeze=1)
        qc_temp = f('quality_index',squeeze=1)
        #y[icount]=xfile[0:4]
        #m[icount]=xfile[4:6]
        #j[icount]=xfile[6:8]
        #h[icount]=xfile[8:10]
        sst_temp.id="NAR SST _ " + xfile[6:8] + xfile[4:6] + xfile[0:4] + " _ " + xfile[8:10] + "h"
        # From Kelvin to Celsius
        sst_temp=numpy.add(sst_temp,-273.15)
        f.close()
        sst_sat.append(sst_temp)
        qc_sst_sat.append(qc_temp)
        ty=int(xfile[0:4])
        tm=int(xfile[4:6])
        tj=int(xfile[6:8])
        th=int(xfile[8:10])
        tt=cdtime.comptime(ty,tm,tj,th)
        time_obs.append(tt)
        icount=icount+1

    
    # Initial VCS:
    #v = vcs.init()
    #data=[lon_obs,lat_obs]
    #v.plot( data )
    #raw_input()

    os.chdir(WORKDIR)

    #return lon_obs, lat_obs, y, m, j, h, sst_sat, qc_sst_sat
    return lon_obs, lat_obs, time_obs, sst_sat, qc_sst_sat
    #-- Fin de lecture des donnees
    #----------------------------------------------------



