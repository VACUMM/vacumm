#!/bin/sh

#*****************************************************************************
#* MOD : SST_evalperf_fig_evolT.sh
#* BUT : Creation de figures d indices statistiques de performance des previsions
#*       pour l evolution temporelle
#* CRE : 24 janvier 2008
#* AUT : P. Craneguy
#* VER : 1.0
#*
#* RQ : SCRIPT_LAND_DETAIL non effectif au 31/01/2008
#*
#****************************************************************************

# Fonctions de datation
. libdate.bash

#-- Entree parametres
#----------------------------------------------------------------------------

if [ $# -ne "14" ]
then
	echo "	      Pas le bon nombre d entrees"	
        exit 1
fi
FILEOUT=${1}
INNOV=${2}
OBS=${3}
BATHY=${4}
OBS_TYPE=${5}	
MODEL_TYPE=${6}
DATE_END=${7}
DUREE=${8}
SCRIPT_LAND_DETAIL=${9}
DIR=${10}
TYPE_FIG=${11}
AREA=${12}
HMAX=${13}
COVMIN=${14}

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

DA1=$Y_DEB$M_D$J_DEB
DA2=$Y_FIN$M_F$J_FIN

#--SCRIPT D EVALUATION DES PERFORMANCES

echo "  --> Evolution temporelle des indices de performance (Figures)"

DIR_FIG=${DIR}/$DA1$sep$DA2
if [ ! -d $DIR_FIG ] ; then
    mkdir $DIR_FIG
fi

FIG0=$DIR_FIG/${OBS_TYPE}_sst_evolT_${AREA}_couverture
FIG1=$DIR_FIG/${OBS_TYPE}_sst_evolT_${AREA}_obs_${MODEL_TYPE}_mean
FIG2=$DIR_FIG/${OBS_TYPE}_sst_evolT_${AREA}_innov_${MODEL_TYPE}
FIG3=$DIR_FIG/${OBS_TYPE}_sst_evolT_${AREA}_rms
FIG4=$DIR_FIG/${OBS_TYPE}_sst_evolT_${AREA}_obs_${MODEL_TYPE}_std
FIG5=$DIR_FIG/${OBS_TYPE}_sst_evolT_${AREA}_nrms

SCRIPT_FERRET=evalperf.jnl
	cat << EOF > ${SCRIPT_FERRET}
	Set Memory /Size=1000
	!!!Set Memory /Size=100
	cancel d /all
	cancel region /all	
	
	!!-- Donnees
	!! ####################
	use "$FILEOUT"
	
	!!-- Regions
	!! ####################		
	go "def_regions.jnl"	
	def sym ZONE=$AREA
	set region/@(\$ZONE)/t="$DATE1":"$DATE2"
	
	
	!-- Limitation a HMAX
	!-- ######################################################
	Let/BAD=99999. BATHY_HM   IF ${BATHY} LE ${HMAX} THEN ${BATHY} ELSE 99999.
	Let/BAD=999. OBS_H1     IF ${BATHY} LE ${HMAX} THEN ${OBS}   ELSE 999.
	Let/BAD=999. INNOV_H1   IF ${BATHY} LE ${HMAX} THEN ${INNOV} ELSE 999.
	
	!-- Couverture donnee effective (i.e. hors masque Terre)
	!-- ######################################################
	Let nb_Terre = BATHY_HM[ X=@nbd, Y=@nbd ]	
	Let ng_data  = OBS_H1[ X=@ngd, Y=@ngd ]
	Let nb_data  = OBS_H1[ X=@nbd, Y=@nbd ]
	Let coverage_mxy = 100.* ng_data / ( nb_data + ng_data - nb_Terre )
	Let covmin_line=${COVMIN} * coverage_mxy / coverage_mxy

	!-- Couverture minimale pour evaluation performances
	!-- ######################################################
	Let/BAD=999.  OBS_HM     IF coverage_mxy GE ${COVMIN} THEN OBS_H1   ELSE 999.
	Let/BAD=999.  INNOV_HM   IF coverage_mxy GE ${COVMIN} THEN INNOV_H1 ELSE 999.
		
	!!-- Indices performances
	!! ######################
	
	Let/unit="degC"  obs_mxy    = OBS_HM[ X=@AVE, Y=@AVE ]
	Let prev       = OBS_HM - INNOV_HM
	Let/unit="degC"  prev_mxy   = prev[ X=@AVE, Y=@AVE ]
	 
	!! Stats difference OBS - MODELE
	Let /unit="degC" innov_mxy  = INNOV_HM[ X=@AVE, Y=@AVE ]
	Let /unit="degC" innov_maxy = innov_mxy + INNOV_HM[ X=@VAR, Y=@VAR ]^0.5
	Let /unit="degC" innov_mixy = innov_mxy - INNOV_HM[ X=@VAR, Y=@VAR ]^0.5	
	
	
	!! Ecart RMS OBS - MODELE
	Let  i2 = INNOV_HM^2.
	Let/unit="degC"   RMSe_mxy = i2[ X=@AVE, Y=@AVE ]^0.5

	!! Ecart centered RMS OBS - MODELE
	Let  ic2 = (INNOV_HM - innov_mxy)^2.
	Let/unit="degC"   RMSe_cent_mxy = ic2[ X=@AVE, Y=@AVE ]^0.5

	!! Variance : OBS & prevision
	Let/unit="degC"  obs_std_mxy  = OBS_HM[ X=@VAR, Y=@VAR ]^0.5
	Let/unit="degC"  prev_std_mxy = prev[ X=@VAR, Y=@VAR ]^0.5
	
	!! NRMSe
	Let NRMSe_mxy = RMSe_mxy / obs_std_mxy
	Let NRMSe_cent_mxy = RMSe_cent_mxy / obs_std_mxy


	!! FIGURES
	!! ####################
	
	!! Parametrage
	Define Symbol OUT ${TYPE_FIG} ! PS ou GIF
	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN Set Mode Metafile
	Define Symbol MODEL ${MODEL_TYPE}       
	

	!!------------ Taux de couverture moyenne spatiale
   	!-- Limites de l axe des ordonnees
   	Let ord1 = 0
   	Let ord2 = 100

   	!-- Plot couverture donnee
	Plot /Nolab /Vlim=\`ord1\`:\`ord2\` /Symb=28 /Size=0.1 /Color=10 /Set   coverage_mxy
   	ppl Title "(\$ZONE) : ${OBS_TYPE}  coverage ($DATE1:$DATE2)"
   	ppl Plot
   	!-- Plot couverture minimale
   	!Plot /Nolab /Symb=20 /Size=0.1 /Color=8 /Over covmin_line 
	Plot /Nolab /Over covmin_line 
	
	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG0}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG0}.gif"
	ENDIF	
			
	!!-------- SST moyenne spatiale
   	!-- Limites de l axe des ordonnees
   	Let ord1 = 11
   	Let ord2 = 30

   	!-- Plot obs
   	Plot /Nolab /Vlim=\`ord1\`:\`ord2\` /Symb=28 /Size=0.1 /Color=10 /Set   obs_mxy
   	ppl Title "(\$ZONE) : mean SST [(^oC, $DATE1:$DATE2)"
   	ppl Plot
   	!-- Plot model
   	Plot /Nolab /Symb=20 /Size=0.1 /Color=8 /Over prev_mxy
   	!-- Legende
   	Let xpct = 0.75  ! Abscisse du coin nord-ouest (de 0 a 1)
   	Let ypct = 0.90  ! Ordonnee du coin nord-ouest (de 0 a 1)
   	Let deltay = 15  ! Espace entre les lignes (1/deltay de la hauteur totale)
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct\`, -1, 0, 0.2, "@c010@PM28@P1@SR  ${OBS_TYPE}"
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct-1*(\$ppl\$ylen)/deltay\`, -1, 0, 0.2, "@P8@PM20@P1@SR  (\$model) "
	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG1}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG1}.gif"
	ENDIF	
	
	!!-------- Innovation moyenne, min, max
   	!-- Limites de l axe des ordonnees
   	Let ord1 = -5
   	Let ord2 = 5

   	!-- Plot obs
   	Plot /Nolab /Vlim=\`ord1\`:\`ord2\` /Symb=28 /Size=0.1 /Color=10 /Set   innov_mxy
   	ppl Title "(\$ZONE) : Bias = ${OBS_TYPE} - ${MODEL_TYPE} (^oC, $DATE1:$DATE2)"
   	ppl Plot
   	!-- Plot model
   	Plot /Nolab /Symb=20 /Size=0.1 /Color=8 /Over innov_mixy
	Plot /Nolab /Symb=32 /Size=0.1 /Color=8 /Over innov_maxy
   	!-- Legende
   	Let xpct = 0.75  ! Abscisse du coin nord-ouest (de 0 a 1)
   	Let ypct = 0.90  ! Ordonnee du coin nord-ouest (de 0 a 1)
   	Let deltay = 15  ! Espace entre les lignes (1/deltay de la hauteur totale)
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct\`, -1, 0, 0.2, "@c010@PM28@P1@SR  Moyenne"
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct-1*(\$ppl\$ylen)/deltay\`, -1, 0, 0.2, "@P8@PM20@P1@SR  - STD "
	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct-2*(\$ppl\$ylen)/deltay\`, -1, 0, 0.2, "@P8@PM20@P1@SR  + STD "
	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG2}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG2}.gif"
	ENDIF	
		
	!!----------- RMS moyenne spatiale
   	!-- Limites de l axe des ordonnees
   	Let ord1 = 0
   	Let ord2 = 5

   	!-- Plot obs
   	Plot /Nolab /Vlim=\`ord1\`:\`ord2\` /Symb=28 /Size=0.1 /Color=10 /Set   RMSe_mxy
   	ppl Title "(\$ZONE) : RMSe SST [ ${OBS_TYPE}, ${MODEL_TYPE} ] (^oC, $DATE1:$DATE2)"
   	ppl Plot
	Plot /Nolab /Symb=20 /Size=0.1 /Color=8 /Over RMSe_cent_mxy
	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG3}.ps metafile.plt
	ENDIF
	!-- Legende
   	Let xpct = 0.75  ! Abscisse du coin nord-ouest (de 0 a 1)
   	Let ypct = 0.90  ! Ordonnee du coin nord-ouest (de 0 a 1)
   	Let deltay = 15  ! Espace entre les lignes (1/deltay de la hauteur totale)
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct\`, -1, 0, 0.2, "@c010@PM28@P1@SR  full RMSe"
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct-1*(\$ppl\$ylen)/deltay\`, -1, 0, 0.2, "@P8@PM20@P1@SR  centered RMSe"
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG3}.gif"
	ENDIF	


	!!----------- STD obs et STD previ
   	!-- Limites de l axe des ordonnees
   	Let ord1 = 0
   	Let ord2 = 5

   	!-- Plot obs
	Plot /Nolab /Vlim=\`ord1\`:\`ord2\`  /Color=10 /Set   obs_std_mxy
   	Plot /Nolab /over /Symb=28 /Size=0.1 /Color=10 /Set   obs_std_mxy
   	ppl Title "(\$ZONE): STD SST [ ${OBS_TYPE}, ${MODEL_TYPE} ] ($DATE1:$DATE2)"
   	ppl Plot
 	!-- Plot model
   	Plot /Nolab /Symb=20 /Size=0.1 /Color=8 /Over prev_std_mxy
   	!-- Legende
   	Let xpct = 0.75  ! Abscisse du coin nord-ouest (de 0 a 1)
   	Let ypct = 0.90  ! Ordonnee du coin nord-ouest (de 0 a 1)
   	Let deltay = 15  ! Espace entre les lignes (1/deltay de la hauteur totale)
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct\`, -1, 0, 0.2, "@c010@PM28@P1@SR  STD ${OBS_TYPE}"
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct-1*(\$ppl\$ylen)/deltay\`, -1, 0, 0.2, "@P8@PM20@P1@SR  STD (\$model) "
	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG4}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG4}.gif"
	ENDIF	
	
	!!--------- NRMS moyenne spatiale

   	!-- Limites de l axe des ordonnees
   	Let ord1 = 0
   	Let ord2 = 5

   	!-- Plot obs
   	Plot /Nolab /Vlim=\`ord1\`:\`ord2\` /Symb=28 /Size=0.1 /Color=10 /Set   NRMSe_mxy
   	ppl Title "(\$ZONE) : NRMSe SST [ ${OBS_TYPE}, ${MODEL_TYPE} ] ($DATE1:$DATE2)"
 	ppl Plot
	Plot /Nolab /Symb=20 /Size=0.1 /Color=8 /Over NRMSe_cent_mxy
	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG3}.ps metafile.plt
	ENDIF
	!-- Legende
   	Let xpct = 0.75  ! Abscisse du coin nord-ouest (de 0 a 1)
   	Let ypct = 0.90  ! Ordonnee du coin nord-ouest (de 0 a 1)
   	Let deltay = 15  ! Espace entre les lignes (1/deltay de la hauteur totale)
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct\`, -1, 0, 0.2, "@c010@PM28@P1@SR  full NRMSe"
   	Label/Nouser \`(\$ppl\$xlen)*xpct\`, \`(\$ppl\$ylen)*ypct-1*(\$ppl\$ylen)/deltay\`, -1, 0, 0.2, "@P8@PM20@P1@SR  centered NRMSe"
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG5}.gif"
	ENDIF	
	
EOF

ferret -nojnl -gif -script ${SCRIPT_FERRET}

#-- Nettoyage
#----------------------------------------------------------------------------
#rm ${SCRIPT_FERRET}
if [ $TYPE_FIG == "PS" ] ; then
	rm metafile*.plt*
fi	
#
#-- copy dans la cible pour visu
#----------------------------------------------------------------------------
#\cp $DIR_FIG/*evolT*.gif /home2/sparta/ftp/SOCOM/EVAL_MENOR/.

exit 0
