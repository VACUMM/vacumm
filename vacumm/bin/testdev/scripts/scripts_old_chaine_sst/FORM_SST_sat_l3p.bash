#!/bin/bash
#**********************************************************************************
#* MOD : FORM_SST_sat_l3p.bash
#* OBJ : Recuperation par FTP anonymous des fichiers de SST L3P sur une 
#        periode de $nSSTdayfile jours depuis la date $DATE_FIN puis concatenation de
#*       ces fichiers dans $FILESORT_L3P      
#* CRE : P. Craneguy (Actimar) / Projet PRECOC (ANR)
#* VER : 1.0 (25/01/2007)
#**********************************************************************************

#-- Arguments
#----------------------------------------------------
# -- passage d arguments
if [ $# -ne 5 ]; then
   echo "usage: submit_L3Psst_FORM.bash $1 $2 $3 $4 $5"
   echo "		HDIR=$1 		# Home dir		"
   echo "		DIR_SST_L3P=$2		# Repertoire de stockage des fichiers de SST (avt concatenation)	"
   echo "		nf=$3			# nbre de jours entre aujourd'hui (date courante) et date de fin de     "
   echo "					# recuperation des donnees sat"
   echo "		nSST=$4			# nbre de jours entre date debut recuperation des donnees sat et  "
   echo "					# aujourd'hui (date courante)"
   echo "		FILESORT_L3P=$5		# Fichier L3P"
   exit 1
fi

HDIR=$1
DIR_SST_L3P=$2
nf=$3
nSST=$4
FILESORT_L3P=$5

#-- arguments fixes
. submit_SST_evalperf_env.bash # Ttes les infos necessaires a recuperation des donnees

#-- Get the path for bash environment
. setenv_NETCDF.bash # Numero de version de NETCDF
source /usr/share/modules/init/bash
module load cmkl/recent
setenv MKL_SERIAL YES
module load netcdf-intel-${NETCDF_INT}/${NETCDF_VERS}


#-- repertoire temporaire (suppr en fin de script)
#----------------------------------------------------
cd $HDIR
TMP_DIR=tmpdir
mkdir ${TMP_DIR}

#----------------------------------------------------
echo ""
echo "---------- RECUPERATION FICHIERS L3P	   -"
echo ""
#----------------------------------------------------


#-- recuperation des donnees par FTP anonymous
#----------------------------------------------------
while [ $nf -le $nSST ]
do 
    an=`date.pl -$nf 3`
    mois=`date.pl -$nf 2`
    jour=`date.pl -$nf 1`
       
    RECUP_FTP_SST_sat_l3p.bash ${L3P_TYPE} ${URL_L3P_DATA} ${FRAC_L3P} $an $mois $jour ${COMPRESS_L3P} ${HDIR}/${DIR_SST_L3P} $NTRY 
  
    nf=$(( ${nf} + 1 ))
done

	
#-- decompression des fichiers journaliers
#----------------------------------------------------
echo  `date` " Decompression des fichiers journaliers "

LISTFILE=`ls $HDIR/${DIR_SST_L3P}/*.nc.${COMPRESS_L3P}`
for FILE in $LISTFILE
do
  bzip2 -d -f $FILE
  if [ $? -ne 0 ]
    then
     rm -f $FILE
  fi
done
	
#-- Reduction a zone UTILE
#----------------------------------------------------
echo  `date` " Reduction a zone UTILE + Time unlimited "
cd $HDIR/${DIR_SST_L3P}
LISTFILE=`ls *.nc`
cd $HDIR
for FILE in $LISTFILE
do
   ncks -d lon,$nxmi,$nxma -d lat,$nymi,$nyma  "${DIR_SST_L3P}/$FILE" "$HDIR/${TMP_DIR}/inc.nc" 
   #-- passage en time unlimited
   ncecat -h -O "$HDIR/${TMP_DIR}/inc.nc"  "$HDIR/${TMP_DIR}/out.nc"
   ncpdq -h -O -a time,record "$HDIR/${TMP_DIR}/out.nc" "$HDIR/${TMP_DIR}/out.nc" 
   ncwa -h -O -a record "$HDIR/${TMP_DIR}/out.nc" "$HDIR/${TMP_DIR}/l3psst_$FILE"
   rm -f "$HDIR/${TMP_DIR}/out.nc" "$HDIR/${TMP_DIR}/inc.nc"
done

#----------------------------------------------------
echo ""
echo "---------- CONCATENATION GRILLE REG	    -"
echo ""
#----------------------------------------------------

cd $HDIR
LISTFILE=`ls ./${TMP_DIR}/l3psst_*.nc`
if [ "$LISTFILE" != "" ]
  then
     ./concat_netcdf.csh "$LISTFILE" "${FILESORT_L3P}"
  else
     echo "--> Aucun fichier a concatener"
fi
  


#------------------------------------------------------------
echo ""
echo "---------- SUPPRESSION ${TMP_DIR} et COMPRESSION	    -"
echo ""
#------------------------------------------------------------

cd $HDIR	
\rm -r $HDIR/${TMP_DIR}

if [ "$SST_DAY_SAVE" == "OUI" ]
    then 
	LISTFILE=`ls ${DIR_SST_L3P}/*.nc`
	if [ "$LISTFILE" != "" ]
	  then
		for FILE in $LISTFILE
		do
  	  		bzip2 $FILE
  	  		if [ $? -ne 0 ]
    	    		then
     				rm -f $FILE
  	  		fi
		done
	fi	
else
   \rm -r ${DIR_SST_L3P}
fi

exit 0
