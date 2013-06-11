#!/bin/sh

#*****************************************************************************
#* MOD : COLLOC_NAR_model.sh
#* BUT : Creation d'un fichiers de donnees NAR et d'innovation (SST_NAR-H.X) 
#        avec H.X=previ modele collocalises dans l'espace (grille OBS)
#* CRE : 04 mars 2008
#* AUT : P. Craneguy
#* VER : 1.0
#****************************************************************************

# Fonctions de datation
. libdate.bash

#-- Entree parametres
#----------------------------------------------------------------------------
if [ $# -ne "21" ]
then
        echo "usage : "
	echo "        $1  FILE_PREV   [fichier de prevision modele] "
	echo "	      $2  VAR_PREV    [nom variable modele] "
	echo "        $3  MVAL_PREV   [missing val de la var modele]"
	echo "        $4  VMIN_PREV   [valeur minimale non acceptee pour la var modele]"
	echo "        $5  VAR_BATHY   [nom variable bathymetrie] "
	echo "        $6  MVAL_ BATHY [missing val de la var bathymetrie]"
	echo "        $7  FILE_OBS    [fichier des donnees] "
	echo "        $8  VAR_OBS     [nom variable obs] "
	echo "        $9  UNIT_OBS    [unite obs KELVIN/CELSIUS] "
	echo "	      ${10} VAR_LEV       [nom variable confidence level] "
	echo "	      ${11} LEV_MIN       [niveau de confiance minimal]"
	echo "	      ${12} FILEOUT       [fichier de sortie]"
	echo "	      ${13} VAR_OUT_PREV  [nom variable sortie OBS - PREVI collocalisee]"
	echo "	      ${14} VAR_OUT_OBS   [nom variable sortie OBS collocalisee]"
	echo "	      ${15} VAR_OUT_BATHY [nom variable sortie bathymetrie] "
	echo "	      ${16} MVAL_OUT      [missing val  des var sorties] " 
	echo "	      ${17} DATE_END	  [Date de fin enregistrement] " 
	echo "	      ${18} DUREE	  [Duree enregistrement en jours] " 
	echo "	      ${19} HEURE	  [Heure de debut d enregistrement] " 
	echo "	      ${20} D_HEURE	  [Pas de temps d enregsitrement en heure] " 
	echo "	      ${21} MODEL_TYPE	  [Type de modele : MARS ou MFS ] " 
        exit 1
fi

FILE_PREV=$1
VAR_PREV=$2
MVAL_PREV=$3
VMIN_PREV=$4
VAR_BATHY=$5
MVAL_BATHY=$6
FILE_OBS=$7
VAR_OBS=$8
UNIT_OBS=${9}
VAR_LEV=${10}
LEV_MIN=${11}
FILEOUT=${12}
VAR_OUT_PREV=${13}
VAR_OUT_OBS=${14}
VAR_OUT_BATHY=${15}
MVAL_OUT=${16}
DATE_END=${17}
DUREE=${18}
HEURE=${19}
D_HEURE=${20}
MODEL_TYPE=${21}

#-- Date de debut de diagnostic
#----------------------------------------------------------------------------
nf=0 
an=`date.pl -$nf 3`
mois=`date.pl -$nf 2`
jour=`date.pl -$nf 1`
JJ=$(date2julien $an $mois $jour) # today julian
Y0=${DATE_END:0:4}
M0=${DATE_END:5:2}
D0=${DATE_END:8:2}
J0=$(date2julien $Y0 $M0 $D0)	  # jour julien date debut
nf=$(($JJ - $J0)) 
nSST=$(($DUREE + $nf -1))
an=`date.pl -$nSST 3`
mois=`date.pl -$nSST 2`
jour=`date.pl -$nSST 1`
DATE_DEB=$an-$mois-$jour

#-- passage en deg C
#----------------------------------------------------------------------------
if [ $UNIT_OBS == "KELVIN" ] 
then
	K2C=273.15
else
	K2C=0
fi

#-- Decomposition des dates
#----------------------------------------------------------------------------
J_DEB=`(echo $DATE_DEB | cut -c 9,10)`
M_D=`(echo $DATE_DEB | cut -c 6,7)`
Y_DEB=`(echo $DATE_DEB | cut -c  1-4)`
J_FIN=`(echo $DATE_END | cut -c 9,10)`
M_F=`(echo $DATE_END | cut -c 6,7)`
Y_FIN=`(echo $DATE_END | cut -c  1-4)`

case $M_D in
    01 ) M_DEB="JAN" ;;
    02 ) M_DEB="FEB" ;;
    03 ) M_DEB="MAR" ;;
    04 ) M_DEB="APR" ;;
    05 ) M_DEB="MAY" ;;
    06 ) M_DEB="JUN" ;;
    07 ) M_DEB="JUL" ;;
    08 ) M_DEB="AUG" ;;
    09 ) M_DEB="SEP" ;;
    10 ) M_DEB="OCT" ;;
    11 ) M_DEB="NOV" ;;
    12 ) M_DEB="DEC" ;;                    
esac

case $M_F in
    01 ) M_FIN="JAN" ;;
    02 ) M_FIN="FEB" ;;
    03 ) M_FIN="MAR" ;;
    04 ) M_FIN="APR" ;;
    05 ) M_FIN="MAY" ;;
    06 ) M_FIN="JUN" ;;
    07 ) M_FIN="JUL" ;;
    08 ) M_FIN="AUG" ;;
    09 ) M_FIN="SEP" ;;
    10 ) M_FIN="OCT" ;;
    11 ) M_FIN="NOV" ;;
    12 ) M_FIN="DEC" ;;                    
esac

sep="-"
DATE1=$J_DEB$sep$M_DEB$sep$Y_DEB
DATE2=$J_FIN$sep$M_FIN$sep$Y_FIN


#-- Formation des fichiers Observation et prevision modele collocalises
#---------------------------------------------------------------------------
#   * 1. Les observations sont :
#		a. selectionnees selon leur niveau de confiance
#		b. interpolees sur la grille 
#			temporelle $DATE_DEB:$DATE_END:$HEURE
#   * 2. Les previsions modele sont interpolees 
#	- en temps : sur la grille temporelle $DATE_DEB:$DATE_END:$HEURE
#	- en espace : sur la grille horizontale des obs
#   * 3. Seuls les pixels renseignes en OBS et en PREVI MODELE sont conserves
#   dans les fichiers de sortie
#
#----------------------------------------------------------------------------

echo "-->  Formation du fichier de données collocalisees "
echo "  	$FILEOUT "

MVAL_INF=999999.

if [ $MODEL_TYPE == "MARS" ] 
then
# La bathy est projette sur la grille des obs
SCRIPT_FERRET=gen_and_save_HX_MARS.jnl
	cat << EOF > ${SCRIPT_FERRET}
	Set Memory /Size=1000
	!!!Set Memory /Size=100
	
	cancel d /all
        !-- Lecture des fichiers
        !   ####################
	Use "$FILE_OBS"       ! Fichier observations
	Use "$FILE_PREV"      ! Fichier previsions modele

	!-- Selection obs selon niveau de confiance et prev selon mval et valeur min
	!   ##########################################################################################
        Let /Bad=$MVAL_OUT obs0   = If ( $VAR_LEV[d=1] ge ${LEV_MIN} ) then ${VAR_OBS}[d=1] else $MVAL_OUT
	Let /Bad=$MVAL_OUT prev0  = If ( ${VAR_PREV}[ d=2 ] ne $MVAL_PREV AND ${VAR_PREV}[ d=2 ] gt ${VMIN_PREV} ) then ${VAR_PREV}[ d=2 ] else $MVAL_OUT
	
	!Let  bath0 = If ( ${VAR_BATHY}[ d=2 ] ne $MVAL_BATHY ) then ${VAR_BATHY}[ d=2 ] else (-1)*$MVAL_OUT
	Let/BAD=$MVAL_INF  bath0 = If ( ${VAR_BATHY}[ d=2 ] ne $MVAL_BATHY ) then ${VAR_BATHY}[ d=2 ] else $MVAL_INF
	
	!-- Definition nouvelles grilles
	!   ############################
	Define Axis /X xnew = x[ d=1, gx=$VAR_OBS ] ! grille x,y = celle des obs
	Define Axis /Y ynew = y[ d=1, gy=$VAR_OBS ]
	Define Grid /X=xnew /Y=ynew gridnew0
	
	Define Axis /T="$DATE1 $HEURE:00:00":"$DATE2 $HEURE:00:00":"${D_HEURE}" tnew		
	Define Grid /X=xnew /Y=ynew /T=tnew gridnew1
	
	!-- Projection sur grille commune
	!   #############################
	
	Let/Title="${VAR_OUT_BATHY}"  ${VAR_OUT_BATHY} = bath0[ g=gridnew0 ]
        Let obs1  = obs0[  g=gridnew1 ]
	Let prev1 = prev0[ g=gridnew1 ]
	
	!-- Conserve les points ou obs ET modele ont des donnees valides
	!   ############################################################
	Let /Bad=$MVAL_OUT  obs2   = IF ( prev1  ne $MVAL_OUT ) then obs1  else $MVAL_OUT
	Let /Bad=$MVAL_OUT  prev2  = IF ( obs2   ne $MVAL_OUT ) then prev1 else $MVAL_OUT
	
        Let /Bad=$MVAL_OUT /Title="${VAR_OUT_OBS}"  /unit="degC" ${VAR_OUT_OBS}  = obs2 - $K2C
	Let /Bad=$MVAL_OUT /Title="${VAR_OUT_PREV}" /unit="degC" ${VAR_OUT_PREV} = ${VAR_OUT_OBS} - prev2

	!-- Sauvegarde ds fichier
	!   #####################
	Save /Clobber /File="$FILEOUT" ${VAR_OUT_OBS}
	Save /Append  /File="$FILEOUT" ${VAR_OUT_PREV}
	Save /Append  /File="$FILEOUT" ${VAR_OUT_BATHY}
EOF

ferret -nojnl -script ${SCRIPT_FERRET}

else
# Gestion des previsions MFS :
#	La bathy n est pas projette sur la grille des obs
#	La grille curviligne de MFS est repassee en grille rectangulaire
SCRIPT_FERRET=gen_and_save_HX_MFS.jnl
	cat << EOF > ${SCRIPT_FERRET}
	Set Memory /Size=10000
	
	cancel d /all
        !-- Lecture des fichiers
        !   ####################
	Use "$FILE_OBS"       ! Fichier observations
	Use "$FILE_PREV"      ! Fichier previsions modele

	!-- Selection obs selon niveau de confiance et prev selon mval et valeur min
	!   ##########################################################################################
        Let /Bad=$MVAL_OUT obs0   = If ( $VAR_LEV[d=1] ge ${LEV_MIN} ) then ${VAR_OBS}[d=1] else $MVAL_OUT
	Let /Bad=$MVAL_OUT prev0  = If ( ${VAR_PREV}[ d=2 ] ne $MVAL_PREV AND ${VAR_PREV}[ d=2 ] gt ${VMIN_PREV} ) then ${VAR_PREV}[ d=2 ] else $MVAL_OUT
	
	!-- Definition nouvelles grilles
	!   ############################
	Define Axis /X xnew = x[ d=1, gx=$VAR_OBS ] ! grille x,y = celle des obs
	Define Axis /Y ynew = y[ d=1, gy=$VAR_OBS ]	
	Define Axis /T="$DATE1 $HEURE:00:00":"$DATE2 $HEURE:00:00":"${D_HEURE}" tnew		
	Define Grid /X=xnew /Y=ynew /T=tnew gridnew1

	!-- Grille curviligne en grille rectangulaire
	!   ##############################################
	! For convenience define variables with the input grid 
	let lonin = NAV_LON[d=2]
	let latin = NAV_LAT[d=2]
	! Define output grid and a variable on the output grid 
	let lonlatout = y[gy=ynew] + x[gx=xnew]
	! Compute the mapping to the rectangular grid
	Let mapp = curv_to_rect_map ( lonin,latin,lonlatout,0.2 )
	
	! Apply the mapping to the data fields 
	Let prev01 = curv_to_rect(prev0[k=1], mapp)
		
	!-- Projection sur grille commune
	!   #############################
	
        Let obs1  = obs0[  g=gridnew1 ]
	Let prev1 = prev01[ g=gridnew1 ]
	
	!-- Conserve les points ou obs ET modele ont des donnees valides
	!   ############################################################
	Let /Bad=$MVAL_OUT  obs2   = IF ( prev1  ne $MVAL_OUT ) then obs1  else $MVAL_OUT
	Let /Bad=$MVAL_OUT  prev2  = IF ( obs2   ne $MVAL_OUT ) then prev1 else $MVAL_OUT
	
        Let /Bad=$MVAL_OUT /Title="${VAR_OUT_OBS}"  /unit="degC" ${VAR_OUT_OBS}  = obs2 - $K2C
	Let /Bad=$MVAL_OUT /Title="${VAR_OUT_PREV}" /unit="degC" ${VAR_OUT_PREV} = ${VAR_OUT_OBS} - prev2

	!-- Sauvegarde ds fichier
	!   #####################
	Save /Clobber /File="$FILEOUT" ${VAR_OUT_OBS}
	Save /Append  /File="$FILEOUT" ${VAR_OUT_PREV}
EOF

#-- lanct d un job car le passage du repere curviligne en repere rectangulaire
#   est couteux en temps calcul
qsub -v SCRIPT_FERRET=${SCRIPT_FERRET} batch_eval

fi


#-- Nettoyage
#----------------------------------------------------------------------------
# rm ${SCRIPT_FERRET}
	
exit 0
