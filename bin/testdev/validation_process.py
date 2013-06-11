#! /usr/bin/env python
# ________________________________________________________________
#
# 
# G. Charria et S. Theetten (11/2009)
# ________________________________________________________________
#
# Base sur la chaine d'evaluation des performances modele en SST
# developpe par P. Craneguy (Actimar) - Projet PRECOC (ANR)
# dans la version 1.0 sur 25/01/2007 
# ________________________________________________________________
import os,sys

print '________________________________________________________________'
print ' Validation de la SST '
print '________________________________________________________________'

# -- Actions et Environnement de SST (Utilisation de la SST NAR)
# ----------------------------------------------------------------
print ' Definition des parametres d environnement '
print '----------------------------------------------------------------'



SCRIPT_DIR=os.getcwd()
os.chdir('../../')
BASENAME=os.getcwd()
WORKDIR = os.path.join(BASENAME,'work')            # repertoire de travail
DIR_SST=os.path.join(WORKDIR,'SST')	# repertoire de stockage des donnees (modele+donnees) concatenees
DIR_FIG=os.path.join(WORKDIR,'FIG')	# racine stockage des figures
DIR_SST_NAR=os.path.join(WORKDIR,'SST_NAR')	# repertoire de stockage des donnees SST NAR
DIR_REGRID=os.path.join(WORKDIR,'REGRID') # repertoire pour l interpolation

# Ajout des repertoire contenant les modules dans le sys.path
sys.path.append(os.path.join(BASENAME,'modules/m0_reading') )
import r_sst_nar
sys.path.append(os.path.join(BASENAME,'modules/m3_visu') )
import v_sst
sys.path.append(os.path.join(BASENAME,'modules/m2_stat') )
import s_basic_stat



lance_form_data_sat='False'		# True: rappatriement des donnees SST (NAR,L3P,)
					# False: Les donnes sont deja dans DIR_SST_NAR


#ZONE_NAR='mocc' # Mediterranee Occidentale a relier a la congif MENOR
ZONE_NAR='gasc' # Gascogne a relier a la config MANGA
#ZONE_NAR='mnord' # Mer du Nord a relier a la config MANGA


andeb='2009'
mdeb='10'  # anciennement 10
jdeb='25'
hdeb='00'

anfin='2009'
mfin='11'  # anciennement 11
jfin='10'
hfin='23'


# -----------------------------------------------------------------
# - Creation des repertoires pour la SST et les figures
os.chdir(WORKDIR)  # on se place dans le repertoire de travail
control_current=os.getcwd()
print "Workdir: %(control_current)s"%vars()

#sys.exit(0)

if os.path.isdir(DIR_SST)==False:
  os.mkdir('SST')

if os.path.isdir(DIR_FIG)==False:
  os.mkdir('FIG')

print "--------------------------------------------------------------------------"
print "--------------------------------------------------------------------------"

# -- Module 0 - Obs: Rapatriement des (ou lien vers) les observations
# corespondantes a la periode consideree
# ----------------------------------------------------------------
print ' -- Lecture des observations satellites --    '
print '----------------------------------------------------------------'
# -- SST satellite (NAR) --
if os.path.isdir(DIR_SST_NAR)==False:
  os.mkdir('SST_NAR')

os.chdir(SCRIPT_DIR)
control_current=os.getcwd()
print "Script_dir: %(control_current)s"%vars()

if lance_form_data_sat=='True':
  r_sst_nar.rappatrie(WORKDIR,DIR_SST_NAR,andeb,mdeb,jdeb,hdeb,anfin,mfin,jfin,hfin,ZONE_NAR)
else:
  print "Les donnees SST NAT doivent se trouver dans le repertoire de travail: %(DIR_SST_NAR)s"%vars()

# -- Lecture des donnees de SST
#lon_obs, lat_obs, y_obs, m_obs, j_obs, h_obs, sst_nar, qc_sst_sat = r_sst_nar.read(WORKDIR,DIR_SST_NAR,andeb,mdeb,jdeb,hdeb,anfin,mfin,jfin,hfin,ZONE_NAR)
lon_obs, lat_obs, time_obs, sst_nar, qc_sst_sat = r_sst_nar.read(WORKDIR,DIR_SST_NAR,andeb,mdeb,jdeb,hdeb,anfin,mfin,jfin,hfin,ZONE_NAR)

# -- Module 0 - Model: Rapatriement des (ou lien vers) resultats
# de simulations correspondates a la periode consideree
# ----------------------------------------------------------------


# -- Module 1: Colocalisation du modele et des observations
# ----------------------------------------------------------------
#sst_nar_int2D = c_coloc_obs_mod.lonlatobs_to_lonlatmodel(WORKDIR,DIR_REGRID,lon_obs,lat_obs,lon_model,lat_model,sst_nar)



# -- Module 2: Calcul des statistiques de validation
# ----------------------------------------------------------------

# ---- QUALITE DE LA MESURE ----
# Calcul de la couverture NAR en % - Evolution temporelle
#narcoverage_time = s_stat_sst_nar.covge_time(lon_obs,lat_obs,sst_nar)

# Calcul de la couverture NAR moyenne en % - Carte
#narcoverage_spatial = s_stat_sst_nar.covge_spat(lon_obs,lat_obs,sst_nar)

# ---- MOYENNES ----
# Calcul de la SST moyenne - Evolution temporelle
obs_sstfcttime = s_basic_stat.spatial_average(lon_obs,lat_obs,sst_nar)
#model_sstfcttime = s_basic_stat.spatial_average(lon_model,lat_model,sst_model)

# ---- ECART-TYPE ----
# Calcul de l'ecart type moyen - Evolution temporelle
obs_stdfcttime = s_basic_stat.spatial_std(lon_obs,lat_obs,sst_nar)
#model_stdfcttime = s_basic_stat.spatial_std(lon_model,lat_model,sst_model)

# Calcul de l'ecart type moyen - Carte
#obs_stdmap = s_basic_stat.map_std(lon_obs,lat_obs,sst_nar)
#model_stdmap = s_basic_stat.map_std(lon_model,lat_model,sst_model)

# ---- BIAIS ----
# Calcul du biais moyen - Evolution temporelle
#bias_time = s_basic_stat.bias_time(lon_obs,lat_obs,sst_nar,lon_model,lat_model,sst_model)

# Calcul du biais moyen - Carte
#bias_spatial = s_basic_stat.bias_spatial(lon_obs,lat_obs,sst_nar,lon_model,lat_model,sst_model)

# Calcul du biais maximum - Carte
#bias_max = s_basic_stat.bias_spatial_max(lon_obs,lat_obs,sst_nar,lon_model,lat_model,sst_model)

# Calcul du biais minimum - Carte
#bias_min = s_basic_stat.bias_spatial_min(lon_obs,lat_obs,sst_nar,lon_model,lat_model,sst_model)

# Calcul du biais systematique - Carte
#syst_bias = s_basic_stat.bias_syst(lon_obs,lat_obs,sst_nar,lon_model,lat_model,sst_model)

# ---- RMS ----
# RMS difference (ou error) - Carte
#rmsdiff = s_basic_stat.rmsdiff_spatial(lon_obs,lat_obs,sst_nar,lon_model,lat_model,sst_model)

# RMS centree - Carte
#rmsnorm = s_basic_stat.rmsnorm(lon_obs,lat_obs,sst_nar,lon_model,lat_model,sst_model)

# RMS difference - Evolution temporelle
rmsdiff_obs = s_basic_stat.rmsdiff_time(lon_obs,lat_obs,sst_nar)
#rmsdiff_model = s_basic_stat.rmsdiff(lon_model,lat_model,sst_model)




# -- Module 3: Visualisation
# ----------------------------------------------------------------
# - 1 - Visualisation des cartes de SST: Observations et Modele
# Exemple de trace pour le premier pas de temps ...
v_sst.map(lon_obs,lat_obs,sst_nar[0],DIR_FIG)

# A developper ...
# - Cartes de differences a chaque pas de temps
#v_sst.mapdiff(lon_obs,lat_obs,sst_nar,lon_model,lat_model,sst_model,DIR_FIG)




# -- Fin de la validation
# ________________________________________________________________
