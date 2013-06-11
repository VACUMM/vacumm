#!/bin/bash
#**********************************************************************************
#* MOD : FORM_SST_mfs.bash
#* OBJ : Formation des fichiers journaliers de SST (modele MFS) sur une
#        periode de nSSTdayfile jours a partir de DATE_FIN
#* EXE : submit_SSTmfs_FORM_MAIN.bash $HDIR $DIR_COMP $FOUT_MOD $nSSTdayfile
#*       sur CAPARMOR
#* CRE : P. Craneguy (Actimar)
#* VER : 1.0 (10/01/2007)
#
#* RQ  : uniquement les champs analyses : pour generaliser a chp forecast
#		- supprimer test existant sur champ analyse
#		- adapter recuperation des champs a cas forecast
#
#**********************************************************************************

echo "--------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------"
echo "-  FORMATION DES FICHIERS DE SST PREVUE POUR EVALUATION 	       	       -"
echo "-                 DES PERF DE MFS			               -"
echo "--------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------"

# Fonctions de datation
. libdate.bash

#-- Arguments
#----------------------------------------------------
# -- passage d arguments
if [ $# -ne 5 ]; then
   echo "usage: submit_SSTmfs_FORM_MAIN.bash $1 $2 $3 $4 $5"
   echo "		HDIR=$1 		# Home dir		"
   echo "		DATE_FIN=$2		# Date de fin de recuperation des donnees"
   echo "		nSSTdayfile=$3		# Durée de recuperation des donnees [jours]"
   echo "		FIELD_TYPE=$4			# Type analyses or forecast"
   echo "		FILESORT=$5		# Fichier de SST concatenee"
   exit 1
fi
HDIR=$1
DATE_FIN=$2
nSSTdayfile=$3
FIELD_TYPE=$4
FILESORT=$5

#-- arguments fixes
. submit_SST_evalperf_env.bash # Ttes les infos necessaires a recuperation des donnees

#-- uniquement champs analyses
if [ "${FIELD_TYPE}" != "${MFS_TYPE1}" ] 
then
   echo "Unknown MFS_TYPE ${MFS_TYPE1}"
   exit 1
fi	   

#-- repertoire temporaire (suppr en fin de script)
#----------------------------------------------------
cd $HDIR
TMP_DIR=tmpdir
mkdir ${TMP_DIR} ${TMP_DIR}_nv
	 
#-- Niveau vertical de la SST PREV (pour NCO)
NZ_SST_PREV_M1=$(( ${MFS_NZ_SST} - 1 ))

#-- Get the path for bash environment (pour NCO)
. setenv_NETCDF.bash # Numero de version de NETCDF
source /usr/share/modules/init/bash
module load cmkl/recent
module load netcdf-intel-${NETCDF_INT}/${NETCDF_VERS}  

	
echo "--------------------------------------------------------------"
echo "- 	 RECUPERATION DES PREVISIONS EN MEDITERRANEE       -"
echo "--------------------------------------------------------------"

if [ -f ${FILESORT} ] ; then
     rm -f ${FILESORT}
fi

#-- periode de recuperation des donnees
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
nf=$(($JJ - $J0)) 		  # nbre de jours entre date fin enregistrement et aujourd'hui
nSST=$(($nSSTdayfile + nf))	  # nbre de jours entre date debut enregistrement et 
				  # aujourd'hui : On recupere un jour de plus que les donnees
				  # sat (car on projette sur les donnees sat) mais
				  # ce dernier jour seule la derniere prevision est recuperee
				  
echo ""
echo "--->   	RECUP SUR CORIOLIS : ${MFS_RESOPE} -"
echo ""

n_manq=0 # compteur de fichiers manquants
while [ $nf -le $nSST ]
do 
    an=`date.pl -$nf 3`
    mois=`date.pl -$nf 2`
    jour=`date.pl -$nf 1`
    Y0_sh=${an:2:2}
    cd ${MFS_RESOPE}
    FILE_MOD=${Y0_sh}${mois}${jour}
    LISTFILE=`ls ${FILE_MOD}*${MFS_FIN}.nc` 
       
    if [ $? == 0 ]
       then    
        echo "	Lien symbolique des previsions du $jour/$mois/$an"
	for FILE in $LISTFILE
	 do
	  ln -s -f ${MFS_RESOPE}/$FILE ${HDIR}/${TMP_DIR}/$FILE
  	done
    else
        echo "	Fichiers du $jour/$mois/$an INEXISTANTS sur CORIOLIS"
    fi 
    cd $HDIR
    
    nf=$(( ${nf} + 1 ))
done


echo "--------------------------------------------------------------"
echo "-		    FORMATION DES PREVISIONS SST 	           -"
echo "--------------------------------------------------------------"

ARP_LONGNAME="time"
ARP_TORIG="1980-01-01 00:00:00"
Y0=${ARP_TORIG:0:4}
M0=${ARP_TORIG:5:2}
D0=${ARP_TORIG:8:2}
J0=$(date2julien $Y0 $M0 $D0)
ARP_UNITS="seconds since $ARP_TORIG"
J_2_S=86400

echo ""
echo "---> EXTRACTION TEMPERATURE SURFACE -"
echo ""
cd $HDIR/${TMP_DIR}
LISTFILE=`ls *.nc`
cd $HDIR
for FILE in $LISTFILE
do
   ncks -v ${MFS_VAR_TEMP},time_counter,nav_lon,nav_lat -d deptht,${NZ_SST_PREV_M1},${NZ_SST_PREV_M1} "$HDIR/${TMP_DIR}/${FILE}" "in.nc" 
   YY=20${FILE:0:2} # on recupere date des donnees a partir du nom du fichier
   MM=${FILE:2:2}
   DD=${FILE:4:2}
   echo " 	Fichier du 	" $DD / $MM / $YY  
   JJ=$(date2julien $YY $MM $DD)
   DT=$(( ( ${JJ} - ${J0} ) * ${J_2_S} ))
   ncap -s "time_counter=time_counter+$DT" in.nc out.nc
   mv -f out.nc in.nc
   ncatted -a long_name,time_counter,o,c,"$ARP_LONGNAME" in.nc
   ncatted -a units,time_counter,o,c,"$ARP_UNITS" in.nc
   mv -f in.nc $HDIR/${TMP_DIR}_nv/mfs_${FILE}
done

echo ""
echo "---> CONCATENATION DES PREVISIONS DE SST dans ${FILESORT} -"
echo ""

cd $HDIR
LISTFILE=`ls ./${TMP_DIR}_nv/mfs_*.nc`

if [ "$LISTFILE" != "" ]
  then
     ./concat_netcdf.csh "$LISTFILE" "${FILESORT}"
  else
     echo "--> Aucun fichier a concatener"
fi



echo ""
echo "--->  SUPPRESSION DES REPERTOIRES TEMPORAIRES 	   	   -"
echo ""
\rm -r ${HDIR}/${TMP_DIR}
\rm -r ${HDIR}/${TMP_DIR}_nv

exit 0
