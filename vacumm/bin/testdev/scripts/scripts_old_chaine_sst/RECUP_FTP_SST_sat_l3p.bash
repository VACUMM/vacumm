#!/bin/bash
#*****************************************************************************************************************
#* MOD : RECUP_FTP_SST_sat_l3p.bash
#* OBJ : Recuperation par FTP d un fichier journalier de SST L3P (soit CMS, soit MERSEA final)
#* CRE : P. Craneguy (Actimar) / Projet PRECOC (ANR)
#* VER : 1.0 (25/01/2007)
#*****************************************************************************************************************
# Fonctions de datation
. libdate.bash

#-- entrees 
if [ $# -ne 9 ]; then
   echo "usage: submit_L3Psst_RECUP_FTP.bash $1 $2 $3 $4 $5 $6 $7 $8 $9	"
   echo "		L3P_TYPE=$1 	# CMS ou MERSEA		"
   echo "		URL_SITE=$2	# Adresse site URL	"
   echo "		EXT=$3		# Nom du fichier : partie invariante"
   echo "		YYYY=$4		# Nom du fichier : année"
   echo "		MM=$5		# Nom du fichier : mois"
   echo "		DD=$6		# Nom du fichier : jour"
   echo "		COMPRESS=$7     # Nom du fichier : extension compression"
   echo "		DIR=$8		# Racine du repertoire de copie du fichier"
   echo "		NTRY=$9		# Nombre d essais de recuperation par FTP"   
   exit 1
fi
L3P_TYPE=$1 	# CMS ou MED
URL_SITE=$2	# Adresse site URL
EXT=$3		# Nom du fichier : partie invariante	
YYYY=$4		# Nom du fichier : année
MM=$5		# Nom du fichier : mois
DD=$6		# Nom du fichier : jour
COMPRESS=$7     # Nom du fichier : extension compression
DIR=$8		# Racine du repertoire de copie du fichier
NTRY=$9		# Nombre d essais de recuperation par FTP


#-- URL site et Nom du fichier a recuperer
url_dir_data=ftp://${URL_SITE}

if   [ "$L3P_TYPE" ==  "CMS" ]
   then
	url_dir=${url_dir_data}/${YYYY}/${MM}
	url_file_def=${YYYY}${MM}${DD}${EXT}.${COMPRESS}
elif [ "$L3P_TYPE" ==  "MERSEA" ]
   then
        J0=$(date2julien $YYYY 01  01)
	JJ=$(date2julien $YYYY $MM $DD)
	JF=$(($JJ - $J0 + 1))
	if [ $JF -lt 100 ]; then
   		DF=0$JF
	fi
	url_dir=${url_dir_data}/${YYYY}/${DF}
	echo $url_dir
	url_file_def=${YYYY}${MM}${DD}-${EXT}.${COMPRESS}
fi

outfile=$DIR/${url_file_def}

#-- Extraction FTP du fichier si necessaire
if  [ ! -f ${outfile} ]
    then
	echo `date` "Downloading L3P file " ${url_file_def}
	wget -t ${NTRY} -O ${outfile} ${url_dir}/${url_file_def}
else
	echo `date` "Pas de downloading : ${outfile} existe deja"
fi
exit 0


