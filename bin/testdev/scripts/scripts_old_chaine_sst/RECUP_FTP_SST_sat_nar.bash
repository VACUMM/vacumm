#!/bin/bash
#**********************************************************************************
#* MOD : RECUP_FTP_SST_sat_nar.bash
#* OBJ : Recuperation par FTP d un fichier de SST NAR horaire
#* CRE : P. Craneguy (Actimar) / Projet PRECOC (ANR)
#* VER : 1.0 (25/01/2007)
#**********************************************************************************

#-- entrees 
if [ $# -ne 9 ]; then
   echo "usage: submit_NARsst_RECUP_FTP.bash $1 $2 $3 $4 $5 $6 $7 $8 $9	"
   echo "		URL_SITE=$1	# Adresse site URL	"
   echo "		ZONE=$2		# Adresse : Zone : mocc "
   echo "		YYYY=$3		# Adresse et Nom du fichier : année"
   echo "		MM=$4		# Nom du fichier : mois"
   echo "		DD=$5		# Nom du fichier : jour"
   echo "		HH=$6		# Nom du fichier : heure"
   echo "		EX=$7     	# Nom du fichier : extension 00.nc.gz"
   echo "		DIR=$8		# Racine du repertoire de copie du fichier"
   echo "		NTRY=$9		# Nombre d essais de recuperation par FTP"   
   
   exit 1
fi

URL_SITE=$1	# Adresse site URL
ZONE=$2		# Adresse : Zone : mocc
YYYY=$3		# Adresse et Nom du fichier : année
MM=$4		# Nom du fichier : mois
DD=$5		# Nom du fichier : jour
HH=$6		# Nom du fichier : heure
EX=$7		# Nom du fichier : extension 00.nc.gz
DIR=$8		# Racine du repertoire de copie du fichier
NTRY=$9		# Nombre d essais de recuperation par FTP

#-- URL site et Nom fichier a recuperer
url_dir_nar=ftp://${URL_SITE}
url_dir=${url_dir_nar}/${ZONE}/${YYYY}
url_file_def=${YYYY}${MM}${DD}${HH}${EX}
outfile=$DIR/${url_file_def}

#-- Extraction du fichier
if  [ ! -f ${outfile} ]
    then
    	echo `date` "Downloading NAR ${ZONE} FILE" ${url_file_def}
	wget -t ${NTRY} -O ${outfile} ${url_dir}/${url_file_def}
else
	echo `date` "Pas de downloading : ${url_dir}/${url_file_def} existe deja"
fi

exit 0


