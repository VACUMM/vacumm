#!/bin/bash
#**********************************************************************************
#* MOD : submit_SST_evalperf_par.bash
#* OBJ : Parametrage des actions d evaluation des perf en SST
#*		SST satellite : NAR   & L3P
#*		SST modele    : MENOR & MFS
#* CRE : P. Craneguy (Actimar) - Projet PRECOC (ANR)
#* VER : 1.0 (25/01/2007)
#**********************************************************************************
#----------------------------------------------------------------------------------
#-- 		DEFINITION DES ACTIONS DE L EVALUATION DES PERFORMANCES
#----------------------------------------------------------------------------------
#--- Actions du script
lance_form_data_sat=OUI		# formation des donnees SST (NAR,L3P,)
lance_form_data_prev=OUI	# formation des donnees SST (previ modeles MARS et MFS)
lance_eval_perf=OUI   		# collocalisation des donnees et eval perf

#----------------------------------------------------------------------------------
#-- 		PARAMETRAGE DES ACTIONS SUR LES DONNEES NAR & L3P et MENOR & MFS
#----------------------------------------------------------------------------------

#-- Actions generales sur les donnees
#-------------------------------------------------------------
nSSTdayfile=590          # duree de comparaison
HOUR=00 	    	# comparaisons obs-previ effectuees a $HOUR h
D_HOUR=24 	    	# comparaisons obs-previ effectuees ttes les $D_HOUR h
#nf=2 			# Nombre de jours avant le jour courant
nf=0                  # Nombre de jours avant le jour courant
an=`date.pl -$nf 3`
mois=`date.pl -$nf 2`
jour=`date.pl -$nf 1`
DATE_FIN=$an-$mois-$jour	# Format YYYY-MM-DD	
FIG_FORMAT=GIF 		# Format des figures PS ou GIF
NTRY=1 		 	# Nbre de tentatives de recuperation d un fichier par FTP anonymous

#-- Zones pour evaluation des performances
#------------------------------------------------------------
n_zone=3		# Nombre de zones avec evaluation des performances
ZONE[0]=MENOR		# Nom des zones d evaluation des performances
ZONE[1]=LIGURE
ZONE[2]=GOL
HMAX_TOT=10000		# Immersion max pour evaluation performances
HMAX_SHELF=$HMAX_TOT
			# Taux minimal de couverture de l observation SST
COV_SST_MIN_MEANT=$(( 100*2/$nSSTdayfile + 1 ))	# pour une moyenne temporelle
COV_SST_MIN_EVOLT=25				# pour une moyenne spatiale

echo " --> Taux de couverture : $COV_SST_MIN_MEANT & $COV_SST_MIN_EVOLT"

#-- Actions sur donnees MFS
#-------------------------------------------------------------
MFS_TYPE=analyses	# analyses : champ analyse ou forecast : champ prevu

#-- Actions sur donnees NAR
#-------------------------------------------------------------
ZONE_NAR=mocc		# zone de definition des donnees source
NAR_DAYNIGHT=NIGHT	# NIGHT, DAY ou ALL
LEV_MIN_COLnarprev=3    # Niveau de confiance minimal des donnees NAR exploitees


#-- Actions sur donnees L3P
#  Les donnees L3 peuvent etre choisies parmi 2 types :
#		TYPE 1 : donnees du CMS Lannion 
#		TYPE 2 : donnees MERSEA
#-------------------------------------------------------------
			#-- Choix type de donnees :
L3P_TYPE=MERSEA        	#    CMS pour CMS Lannion
			#    MERSEA pour Mearsea 
		  	# ATTENTION : 
			# Au 30/01/2008 Il n'y avait pas de donnees MERSEA
			# disponible sur le site FTP avant le 14/01/2008 	
