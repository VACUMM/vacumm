#!/bin/bash
#**********************************************************************************
#* MOD : FORM_SST_sat_nar.bash
#* OBJ : Recuperation par FTP anonymous des fichiers de SST NAR sur une 
#        periode de $ndayfile jours depuis la date $DATE_FIN puis concatenation de
#*       ces fichiers dans $FILESORT_NAR     
#* CRE : P. Craneguy (Actimar) / Projet PRECOC (ANR)
#* VER : 1.0 (25/01/2007)
#**********************************************************************************

#-- Arguments
#----------------------------------------------------
# -- passage d arguments
if [ $# -ne 5 ]; then
   echo "usage: submit_NARsst_FORM.bash $1 $2 $3 $4 $5"
   echo "		HDIR=$1 		# Home dir		"
   echo "		DIR_SST_NAR=$2		# Repertoire de stockage des fichiers de SST (avt concatenation)	"
   echo "		nf=$3			# nbre de jours entre aujourd'hui (date courante) et date de fin de     "
   echo "					# recuperation des donnees sat"
   echo "		nSST=$4			# nbre de jours entre date debut recuperation des donnees sat et  "
   echo "					# aujourd'hui (date courante)"
   echo "		FILESORT_NAR=$5		# Fichier NAR"
   exit 1
fi

HDIR=$1
DIR_SST_NAR=$2
nf=$3
nSST=$4
FILESORT_NAR=$5

#-- arguments fixes
source environnement.sh # Ttes les infos necessaires a recuperation des donnees


#-- repertoire temporaire (suppr en fin de script)
#----------------------------------------------------
TMP_DIR=tmpdir
cd $HDIR
mkdir ${TMP_DIR}

#----------------------------------------------------
echo ""
echo "---------- RECUPERATION FICHIERS NAR	   -"
echo ""
#----------------------------------------------------

#-- recuperation de la grille NAR
#----------------------------------------------------
url_dir=ftp://${URL_NAR_GRID}
echo `date` "Downloading NAR ${ZONE_NAR} grid FILE"
wget -O $HDIR/${gridfile}.gz ${url_dir}/${gridfile}.gz
gunzip -f $HDIR/${gridfile}.gz

#-- recuperation des donnees par FTP anonymous
#----------------------------------------------------
while [ $nf -le $nSST ]
do 
    an=`date.pl -$nf 3`
    mois=`date.pl -$nf 2`
    jour=`date.pl -$nf 1`
    
    if [ ${NAR_DAYNIGHT} == NIGHT ] || [ ${NAR_DAYNIGHT} == ALL ] 
    then
      echo " 	--> donnees de nuit"
      H1=${HEURE_NAR[0]}
      H2=${HEURE_NAR[3]}
      RECUP_FTP_SST_sat_nar.bash ${URL_NAR_DATA} ${ZONE_NAR} $an $mois $jour $H1 ${EXT_NAR} ${HDIR}/${DIR_SST_NAR} $NTRY
      RECUP_FTP_SST_sat_nar.bash ${URL_NAR_DATA} ${ZONE_NAR} $an $mois $jour $H2 ${EXT_NAR} ${HDIR}/${DIR_SST_NAR} $NTRY
    fi
    if [ ${NAR_DAYNIGHT} == DAY  ] || [ ${NAR_DAYNIGHT} == ALL ] 
    then
      echo " 	--> donnees de jour"
      H1=${HEURE_NAR[1]}
      H2=${HEURE_NAR[2]}      
      RECUP_FTP_SST_sat_nar.bash ${URL_NAR_DATA} ${ZONE_NAR} $an $mois $jour $H1 ${EXT_NAR} ${HDIR}/${DIR_SST_NAR} $NTRY
      RECUP_FTP_SST_sat_nar.bash ${URL_NAR_DATA} ${ZONE_NAR} $an $mois $jour $H2 ${EXT_NAR} ${HDIR}/${DIR_SST_NAR} $NTRY
    fi        
    nf=$(( ${nf} + 1 ))
done


#----------------------------------------------------
echo ""
echo "---------- MISE SUR GRILLE REGULIERE	   -"
echo ""
#----------------------------------------------------

#-- Decompression des fichiers horaires
#----------------------------------------------------
echo  `date` " Decompression des fichiers horaires "

LISTFILE=`ls $HDIR/${DIR_SST_NAR}/*.nc.${COMPRESS_NAR}`
for FILE in $LISTFILE
do
  gunzip -f $FILE
  if [ $? -ne 0 ]
    then
     rm -f $FILE
  fi
done	
	
#-- fichiers horaires sur grille reguliere
#----------------------------------------------------

cd $HDIR/${DIR_SST_NAR}
LISTFILE=`ls *.nc`
ZONEUTIL="-R latMin=$lami,latMax=$lama,lonMin=$lomi,lonMax=$loma -I latRes=$lare,lonRes=$lore"
for FILE in $LISTFILE
do
   cd $HDIR
   FILEIN="$HDIR/${DIR_SST_NAR}/$FILE"
   FILEOUT="./${TMP_DIR}/narsst_$FILE"
   nar2cf_caparmor -i $FILEIN -g $gridfile -o $FILEOUT $ZONEUTIL
done


#----------------------------------------------------
echo ""
echo "---------- CONCATENATION GRILLE REG	    -"
echo ""
#----------------------------------------------------

cd $HDIR
LISTFILE=`ls ./${TMP_DIR}/narsst_*.nc`
if [ "$LISTFILE" != "" ]
  then
  	./concat_netcdf.csh "$LISTFILE" "${FILESORT_NAR}"
  else
     echo "--> Aucun fichier a concatener"
fi

#------------------------------------------------------------
echo ""
echo "---------- SUPPRESSION ${TMP_DIR} et COMPRESSION	    -"
echo ""
#------------------------------------------------------------

cd $HDIR
\rm -r $HDIR/${TMP_DIR} $HDIR/${gridfile}
if [ "$SST_DAY_SAVE" == "OUI" ]
    then 
	LISTFILE=`ls ${DIR_SST_NAR}/*.nc`	
	if [ "$LISTFILE" != "" ]
	  then
		for FILE in $LISTFILE
		do
  	  		gzip -f $FILE
  	  		if [ $? -ne 0 ]
    	    		then
     				rm -f $FILE
  	  		fi
		done
	fi	
else
   \rm -r ${DIR_SST_NAR}
fi
cd $HDIR		


exit 0
