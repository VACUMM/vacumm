#!/bin/bash
#**********************************************************************************
#* MOD : submit_SST_evalperf_env.bash
#* OBJ : Environnement pour l evaluation des performances en SST du 
#*       modele operationnel MENOR
#*		SST satellite : NAR   & L3P
#*		SST modele    : MENOR & MFS
#* CRE : P. Craneguy (Actimar) - Projet PRECOC (ANR)
#* VER : 1.0 (25/01/2007)
#**********************************************************************************

#-- parametres des actions
. submit_SST_evalperf_par.bash

#----------------------------------------------------------------------------------
#-- 		ACTIONS SUR LES DONNEES NAR & L3P et MENOR & MFS
#----------------------------------------------------------------------------------

#-- Actions generales sur les donnees
#-------------------------------------------------------------
HDIR=$PWD   		# repertoire de travail
DIR_SST=$HDIR/SST		# repertoire de stockage des donnees (sat+prev) concatenees
DIR_FIG=$HDIR/FIG	# racine stockage des figures
SST_DAY_SAVE=NON 	# OUI : on garde les donnees sources de NAR et L3P en compresse
                 	# NON : on les supprime

#-- Actions sur donnees NAR
#-------------------------------------------------------------
FOUT_NAR=SST_NAR.nc
lami=39.5		# grille projection donnees NAR
lama=44.5		# correspond a grille MENOR
lomi=0.0
loma=16.0
lare=.02
lore=.03


#-- Actions sur donnees L3P
#  Les donnees L3 peuvent etre choisies parmi 2 types :
#		TYPE 1 : donnees du CMS Lannion 
#		TYPE 2 : donnees MERSEA
#-------------------------------------------------------------
FOUT_L3P=SST_L3P.nc
			#-- Choix type de donnees :
L3P_TYPE1=CMS     	#    CMS pour CMS Lannion
L3P_TYPE2=MERSEA	#    MERSEA pour Mearsea 

if   [ "$L3P_TYPE" == "${L3P_TYPE1}" ]	#-- Indices de grille pour 
   then					#  correspondance avec grille MENOR 	
	nxmi=2000	 #-- BORNES (ncks) pour MENOR
	nxma=2319	
	nymi=2190
	nyma=2289	
elif [ "$L3P_TYPE" == "${L3P_TYPE2}" ] 
   then
	nxmi=1175	 #-- BORNES (ncks) pour MENOR
	nxma=1974	
	nymi=725
	nyma=974	
else
	echo "Unknown L3P Data Type"
	exit 1
fi

#-- Actions sur les donnees modele
#-------------------------------------------------------------
FOUT_PREV_MARS=SST_prev_mars.nc
FOUT_PREV_MFS=SST_prev_mfs.nc

#-- Actions de collocalisation avec les donnees NAR
#-------------------------------------------------------------
FILEOUT_COLnarprev_NAR=$HDIR/$DIR_SST/NARprev_colloc.nc # modele MARS
FILEOUT_COLnarmfs_NAR=$HDIR/$DIR_SST/NARmfs_colloc.nc   # modele MFS
VAR_OUT_COLnarprev_prev=innov_nar
VAR_OUT_COLnarprev_OBS=sst_nar
VAR_OUT_COLnarprev_bathy=H0_nar

#-- Actions de collocalisation avec les donnees L3P
#-------------------------------------------------------------
#FILEOUT=$HDIR/$DIR_SST/NARprev_colloc.nc
#VAR_OUT_PREV=innov_nar
#VAR_OUT_OBS=sst_nar

MVAL_OUT=999



#----------------------------------------------------------------------------------
#-- 		GENERALITES SUR LES DONNEES NAR & L3P et MENOR & MFS
#----------------------------------------------------------------------------------

#-- Generalites (site URL, type de donnes)
#-------------------------------------------------------------
URL_CERSAT=ftp.ifremer.fr/pub/ifremer/cersat
SST_A=NAR
SST_B=L3P
SCRIPT_LAND_DETAIL=medi_coast

#-- Données NAR
#-------------------------------------------------------------
#------- recuperation des donnees : generalites (site ftp, etc)
URL_NAR_DATA=${URL_CERSAT}/SAFOSI/Products/NARSST/netcdf
URL_NAR_GRID=${URL_CERSAT}/SAFOSI/Products/NARSST/netcdf/grids
gridfile=grid_${ZONE_NAR}.nc
HEURE_NAR=( 02 10 12 20 )
COMPRESS_NAR=gz
EXT_NAR=00.nc.${COMPRESS_NAR}

#-- exploitation donnees pour eval performances
VAR_SST_NAR=temp
UNIT_SST_NAR=degC
VAR_LEV_SST_NAR=conf_level

#-- Données L3P
#-------------------------------------------------------------
#------- recuperation des donnees : generalites (site ftp, etc)
if [ "${L3P_TYPE}" == "${L3P_TYPE1}" ] 
   then
	echo " DONNEES L3P supercollated du CMS"
	URL_L3P_DATA=${URL_CERSAT}/products/gridded/experimental-cms/netcdf/mersea/supercollated
	FRAC_L3P=00-EUR-scol005x005-EUR-v01.nc
elif [ "${L3P_TYPE}" == "${L3P_TYPE2}" ] 
   then
	echo " DONNEES L3P supercollated de MEDSPIRATION"
	URL_L3P_DATA=ftp.ifremer.fr/ifremer/medspiration/data/l3/eurdac/med/odyssea
	FRAC_L3P=IFREMER-L3_GHRSST-SSTsubskin-MERGED_MED002-adjusted-0000-v02-fv01.nc

else
	echo "Unknown L3P Data Type"
	exit 1
fi
COMPRESS_L3P=bz2
VDIM_L3P=( time lat lon )

#------- exploitation donnees pour eval performances

		
#-- Modele MENOR 
#-------------------------------------------------------------
#------- recuperation des donnees :
MARS_CONF1_PREV_ENTETE=champs_V7.42_
# entete fichier resultat : menor (avt dec 2007) menor_ ensuite
VAR_HO_PREV=H0		# variable H0 ds fichier de resultat
VAR_TEMP_PREV=TEMP	# variable temperature ds fichier de resultat
NZ_SST_PREV=30		# niv vertical de la SST
LAST_HOUR_PREV=21	# derniere haure prevue	
#	-- Previsions sur CORIOLIS -----------------------------------------------
#	RESOPE_PREV=/home/coriolis_exp/spool/co01/co0123/co012302/co01230207/f2
        RESOPE_PREV=/home2/creizic/pgarreau/MENOR_2007
#	-- Previsions sur RES      -----------------------------------------------
	RESRUN_RAC=/home13/caparmor/previmer/op/socom/res #--Pierre Garreau
	MARS_CONF1_PREV=MENOR	
        MARS_CONF2_PREV=${MARS_CONF1_PREV}-V6.18
	RESRUN_PREV=${RESRUN_RAC}/${MARS_CONF1_PREV}/${MARS_CONF2_PREV}
	
#------- exploitation donnees pour eval performances
VAR_SST_PREV=temp
MVAL_SST_PREV=99
VMIN_SST_PREV=0 # valeur minimale de VAR_SST_PREV
MVAL_H0_PREV=$(((-1)*$MVAL_SST_PREV))

#-- Modele MFS
#-------------------------------------------------------------
MFS_FIN=_T		# fin du nom de fichier MFS pour salinite et tempreature
MFS_VAR_TEMP=votemper   # meme nom pour variable forecast
MFS_NZ_SST=1		# SST : 1er niveau vertical
MFS_TYPE1=analyses	# champ analyse MFS
MFS_TYPE2=forecast	# champ prevu MFS

#	-- Previsions sur CORIOLIS -----------------------------------------------
        MFS_RESOPE_RAC=/home/coriolis_exp/spool/co01/co0123/co012302/co01230209/archive
	if [ "${MFS_TYPE}" == "${MFS_TYPE1}" ] 
   	then
	   MFS_RESOPE=${MFS_RESOPE_RAC}/${MFS_TYPE1}
	else
	   echo "Unknown MFS_TYPE ${MFS_TYPE1}"
	   exit 1
        fi	   
	
#	-- Previsions sur CORIOLIS
