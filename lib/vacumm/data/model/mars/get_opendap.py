#utilisation de ncread_best_estimate

# exemple avec fichier a l'IDRIS.
#xfile='http://dods.idris.fr/reee512/FORCINGS/ORCA_R2/NEMO/hum_cloud_1m.nc'
#var = ncread_best_estimate("socliohu",xfile)

# ou 

#xfile='http://www.ifremer.fr/thredds3/standard/fileServer/PREVIMER_F1-MARS3D-MANGA/2008/PREVIMER_F1-MARS3D-MANGA_20080108T0000Z.nc'

# ou
# xfile ='http://www.ifremer.fr/thredds3/standard/dodsC/PREVIMER_F1-MARS3D-MANGA_FULL_TIME_SERIE'

#ou
# xfile='http://www.ifremer.fr/thredds3/standard/dodsC/PREVIMER_F1-MARS3D-MANGA/2008/PREVIMER_F1-MARS3D-MANGA_20080108T0100Z.nc'

# ==> erreur due manifestement a cdms.open ....
# en effet :
#f=cdms2.open(xfile)
# ne fonctionne pas...
