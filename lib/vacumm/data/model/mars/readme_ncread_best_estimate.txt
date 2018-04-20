# utilisation de ncread_best_estimate

# exemples d'utilisation :
# lecture d'un seul fichier et chargement de la variable TEMP 
# model1 = ncread_best_estimate('TEMP','PREVIMER_F1-MARS3D-MANGA_20091026T0700Z.nc')

# lecture de plusieurs fichiers avec metacaracteres dans le nom et chargement de la variable TEMP 

# attention, il peut y avoir des erreurs suite a un mauvais tri des fichiers par date.
# par exemple, si les dates ont un format de type 01-01-2010, ca les inversera ....
# model2 = ncread_best_estimate('TEMP','PREVIMER_F1-MARS3D-MANGA_200910??T0700Z.nc')
# il vaut mieux definir un file_pattern avec l'option time

# lecture de plusieurs fichiers entre ctdeb et ctfin et chargement de la variable TEMP 
# model3 = ncread_best_estimate('TEMP',"PREVIMER_F1-MARS3D-MANGA_%Y%m%dT0700Z.nc", (ctdeb, ctfin))
