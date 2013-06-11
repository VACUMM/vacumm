#!/bin/bash
#**********************************************************************************
#* MOD : FORM_SST_sat_main.bash
#* OBJ : Recuperation par FTP anonymous des fichiers de SST NAR et L3P sur une 
#        periode de $nSSTdayfile jours depuis la date $DATE_FIN puis concatenation de
#*       ces fichiers dans $FILESORT_NAR et $FILESORT_L3P      
#* CRE : P. Craneguy (Actimar) - Projet PRECOC (ANR)
#* VER : 1.0 (25/01/2007)
#**********************************************************************************

# Fonctions de datation
. libdate.bash


#-- Arguments
#------------------------------------
# -- passage d arguments
if [ $# -ne 7 ]; then
   echo "usage: submit_SSTsat_FORM_MAIN.bash $1 $2 $3 $4 $5 $6 $7"
   echo "		HDIR=$1 		# Home dir		"
   echo "		DIR_SST_NAR=$2	 	# Repertoire de stockage des fichiers de SST NAR (avt concatenation)"
   echo "		DIR_SST_L3P=$3	 	# Repertoire de stockage des fichiers de SST L3P (avt concatenation)"
   echo "		FILESORT_NAR=$4		# Fichier concatene NAR (avec chemin complet)"
   echo "		FILESORT_L3P=$5		# Fichier concatene L3P (avec chemin complet)"
   echo "		DATE_FIN=$6		# Date de fin de recuperation des donnees"
   echo "		nSSTdayfile=$7		# Durée de recuperation des donnees [jours]"
   exit 1
fi

HDIR=$1		 
DIR_SST_NAR=$2	 
DIR_SST_L3P=$3	 
FILESORT_NAR=$4	 
FILESORT_L3P=$5	 
DATE_FIN=$6	 
nSSTdayfile=$7  

#-- arguments fixes
. submit_SST_evalperf_env.bash # Ttes les infos necessaires a recuperation des donnees

#------------------------------------	

#-- periode de recuperation des donnees satellites
#----------------------------------------------------
nf=0 
an=`date.pl -$nf 3`
mois=`date.pl -$nf 2`
jour=`date.pl -$nf 1`
JJ=$(date2julien $an $mois $jour) # today julian
Y0=${DATE_FIN:0:4}
M0=${DATE_FIN:5:2}
D0=${DATE_FIN:8:2}
J0=$(date2julien $Y0 $M0 $D0)	  # jour julien date debut
nf=$(($JJ - $J0)) 		  # nbre de jours entre aujourd'hui (date courante) et date de fin de
				  # recuperation des donnees sat
nSST=$(($nSSTdayfile + nf -1))	  # nbre de jours entre date debut recuperation des donnees sat et  
				  # aujourd'hui (date courante)
			  		
echo "------------------------------------------------------"
echo "-  	DONNEES  NAR : recup & formation	   -"
echo "------------------------------------------------------"

if [ -f ${FILESORT_NAR} ] ; then
     rm -f ${FILESORT_NAR}
fi

FORM_SST_sat_nar.bash ${HDIR} ${DIR_SST_NAR} ${nf} ${nSST} ${FILESORT_NAR}


##echo "------------------------------------------------------"
##echo "-  	DONNEES  L3P : recup & formation	   -"
##echo "------------------------------------------------------"
##
##if [ -f ${FILESORT_L3P} ] ; then
##     rm -f ${FILESORT_L3P}
##fi	
##
##FORM_SST_sat_l3p.bash ${HDIR} ${DIR_SST_L3P} ${nf} ${nSST} ${FILESORT_L3P}
##
exit 0
		
