#!/bin/bash
#**********************************************************************************
#* MOD : submit_SSTprev_FORM_MAIN.bash
#* OBJ : Formation des fichiers journaliers de SST (modele MENOR) sur une
#        periode de nSSTdayfile jours a partir de DATE_FIN
#* EXE : submit_formSST_sat_main.csh $HDIR $DIR_COMP $FOUT_MOD $nSSTdayfile
#*       sur CAPARMOR
#* CRE : P. Craneguy (Actimar)
#* VER : 1.0 (10/01/2007)
#**********************************************************************************

echo "--------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------"
echo "-  FORMATION DES FICHIERS DE SST PREVUE POUR EVALUATION 	       	       -"
echo "-                 DES PERF DE MARS			               -"
echo "--------------------------------------------------------------------------"
echo "--------------------------------------------------------------------------"

# Fonctions de datation
. libdate.bash

#-- Arguments
#----------------------------------------------------
# -- passage d arguments
if [ $# -ne 4 ]; then
   echo "usage: submit_SSTprev_FORM_MAIN.bash $1 $2 $3 $4"
   echo "		HDIR=$1 		# Home dir		"
   echo "		DATE_FIN=$2		# Date de fin de recuperation des donnees"
   echo "		nSSTdayfile=$3		# Durée de recuperation des donnees [jours]"
   echo "		FILESORT_PREV=$4	# Fichier de SST concatenee"
   exit 1
fi
HDIR=$1
DATE_FIN=$2
nSSTdayfile=$3
FILESORT_PREV=$4

#-- arguments fixes
. submit_SST_evalperf_env.bash # Ttes les infos necessaires a recuperation des donnees


#-- repertoire temporaire (suppr en fin de script)
#----------------------------------------------------
cd $HDIR
TMP_DIR=tmpdir
mkdir ${TMP_DIR} ${TMP_DIR}_nv
	 
#-- Niveau vertical de la SST PREV (pour NCO)
NZ_SST_PREV_M1=$(( ${NZ_SST_PREV} - 1 ))

#-- Get the path for bash environment (pour NCO)
. setenv_NETCDF.bash # Numero de version de NETCDF
source /usr/share/modules/init/bash
module load cmkl/recent
module load netcdf-intel-${NETCDF_INT}/${NETCDF_VERS}  

	
echo "--------------------------------------------------------------"
echo "- 	 RECUPERATION DES PREVISIONS EN MEDITERRANEE      -"
echo "--------------------------------------------------------------"

if [ -f ${FILESORT_PREV} ] ; then
     rm -f ${FILESORT_PREV}
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
echo "--->   	RECUP SUR CORIOLIS : ${RESOPE_PREV} -"
echo ""

n_manq=0 # compteur de fichiers manquants
while [ $nf -le $nSST ]
do 
    an=`date.pl -$nf 3`
    mois=`date.pl -$nf 2`
    jour=`date.pl -$nf 1`
    
    cd ${RESOPE_PREV}
    FILE_MOD=${MARS_CONF1_PREV_ENTETE}$an$mois$jour
 echo "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"
 pwd
 echo $FILE_MOD
 echo "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW"
#    if [ $nf == $nSST ] 
#      then
#       LISTFILE=`ls ${FILE_MOD}*${LAST_HOUR_PREV}*.nc*` 
#    else
       LISTFILE=`ls ${FILE_MOD}*.nc*` 
#    fi 
       
    if [ $? == 0 ]
       then    
        echo "	Lien symbolique des previsions du $jour/$mois/$an"
	for FILE in $LISTFILE
	 do
	  ln -s -f ${RESOPE_PREV}/$FILE ${HDIR}/${TMP_DIR}/$FILE
  	done
    else
        echo "	Fichiers de prevision du $jour/$mois/$an INEXISTANTS sur CORIOLIS"
	echo "  --> ces fichiers seront recherches sur RES_MARS"   
        n_manq=$(( ${n_manq} + 1 )) 
	n_manq_M1=$(( ${n_manq} - 1 ))
	nf_manq[$n_manq_M1]=$nf  # jours manquants sur Coriolis       
    fi 
    cd $HDIR
    
    nf=$(( ${nf} + 1 ))
done


#-- recuperation sur RES_MAR des fichiers manquants sur CORIOLIS
if [ $n_manq -ne 0 ]
 then
 echo "--------------------------------------------------------------------"
 echo "    On a $n_manq Fichiers manquants sur CORIOLIS			   "
 echo "--->    RECUPERATION des Fichiers manquants sur ${RESRUN_PREV}     -"
 echo "--------------------------------------------------------------------"
 n_manq_res=0
  for nf in ${nf_manq[*]} # boucle sur jours manquants
    do
     an=`date.pl -$nf 3`
     mois=`date.pl -$nf 2`
     jour=`date.pl -$nf 1`     
     cd ${RESRUN_PREV}
     FILE_MOD=${MARS_CONF1_PREV_ENTETE}$an$mois$jour

 #    if [ $nf == $nSST ] 
 #     then
 #      LISTFILE=`ls ${FILE_MOD}*${LAST_HOUR_PREV}*.nc*` 
 #    else
       LISTFILE=`ls ${FILE_MOD}*.nc*` 
 #    fi      
     if [ $? == 0 ]
       then    
        echo "	Lien symbolique des previsions du $jour/$mois/$an "
	for FILE in $LISTFILE
	 do
	  ln -s -f ${RESRUN_PREV}/$FILE ${HDIR}/${TMP_DIR}/$FILE
  	done	
     else
        echo "	Fichiers de prevision du $jour/$mois/$an INEXISTANTS sur RES_MARS"
        n_manq_res=$(( ${n_manq_res} + 1 )) 
     fi 
     cd $HDIR
  done
  if [ $n_manq_res -ne 0 ] ; then    
      n_recup=$(( ${n_manq} - ${n_manq_res}))
      echo "    On a recupere $n_recup Fichiers sur RES_MARS" 
      echo "    Il manque toutefois $n_manq_res "   
  else
      echo "    --> Tout a ete recupere"     
  fi 
else
  echo ""
  echo "---> Aucun fichier journalier manquant sur CORIOLIS : pas de recup sur RES_MARS"
  echo ""
fi

echo "--------------------------------------------------------------"
echo "-		    FORMATION DES PREVISIONS SST 	           -"
echo "--------------------------------------------------------------"

echo ""
echo "---> EXTRACTION TEMPERATURE SURFACE -"
echo ""
cd $HDIR/${TMP_DIR}
LISTFILE=`ls *.nc*`
cd $HDIR
for FILE in $LISTFILE
do
   #echo ncks -v ${VAR_HO_PREV},${VAR_TEMP_PREV} -d level,${NZ_SST_PREV_M1},${NZ_SST_PREV_M1} "$HDIR/${TMP_DIR}/${FILE}" "$HDIR/${TMP_DIR}_nv/socom_${FILE}" 
   ncks -v ${VAR_HO_PREV},${VAR_TEMP_PREV} -d level,${NZ_SST_PREV_M1},${NZ_SST_PREV_M1} "$HDIR/${TMP_DIR}/${FILE}" "$HDIR/${TMP_DIR}_nv/socom_${FILE}" 
done

echo ""
echo "---> CONCATENATION DES PREVISIONS DE SST dans ${FILESORT_PREV} -"
echo ""

cd $HDIR
LISTFILE=`ls ./${TMP_DIR}_nv/socom_*.nc*`

if [ "$LISTFILE" != "" ]
  then
     ./concat_netcdf.csh "$LISTFILE" "${FILESORT_PREV}"
  else
     echo "--> Aucun fichier a concatener"
fi

echo ""
echo "--->  SUPPRESSION DES REPERTOIRES TEMPORAIRES 	   	   -"
echo ""
\rm -r ${HDIR}/${TMP_DIR}
\rm -r ${HDIR}/${TMP_DIR}_nv

exit 0
