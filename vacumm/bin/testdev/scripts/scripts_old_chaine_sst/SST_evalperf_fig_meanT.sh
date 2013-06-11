#!/bin/sh

#*****************************************************************************
#* MOD : SST_evalperf_fig.sh
#* BUT : Creation de figures d indices statistiques de performance des previsions
#*       		en moyenne temporelle
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

echo "  --> Distribution spatiale des indices de performance (Figures)"

DIR_FIG=${DIR}/$DA1$sep$DA2
if [ ! -d $DIR_FIG ] ; then
    mkdir $DIR_FIG
fi

FIG0=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_couverture
FIG1=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_obs
FIG2=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_${MODEL_TYPE}
FIG3=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_innov
FIG4=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_rmse	
FIG4bis=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_rmse_centered	
FIG5=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_innov_min
FIG6=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_innov_max
FIG7=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_innov_bias
FIG8=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_std_obs
FIG9=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_std_${MODEL_TYPE}
FIG10=$DIR_FIG/${OBS_TYPE}_sst_meanT_${AREA}_std_obs${MODEL_TYPE}

SCRIPT_FERRET=evalperf.jnl
	cat << EOF > ${SCRIPT_FERRET}
	Set Memory /Size=1000
	!!!Set Memory /Size=100
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
	Let/BAD=99999. BATHY_HM  IF ${BATHY} LE ${HMAX} THEN ${BATHY} ELSE 99999.
	Let/BAD=999.   OBS_H1    IF ${BATHY} LE ${HMAX} THEN ${OBS}   ELSE 999.
	Let/BAD=999.   INNOV_H1  IF ${BATHY} LE ${HMAX} THEN ${INNOV} ELSE 999.
	!Let/BAD=999.   OBS_HM    IF ${BATHY} LE ${HMAX} THEN ${OBS}   ELSE 999.
	!Let/BAD=999.   INNOV_HM  IF ${BATHY} LE ${HMAX} THEN ${INNOV} ELSE 999.
	
	!!-- Masque Terre
	!! ####################	
	let mask = if BATHY_HM NE -999. THEN 0 ELSE 1
	let mask_terre = if mask EQ 1 then 1

	!-- Couverture donnee effective (i.e. hors masque Terre)
	!-- ######################################################
	Let coverage_mt1  = 100.*OBS_H1[T=@ngd] / ( OBS_H1[T=@ngd] + OBS_H1[T=@nbd] )
	!Let coverage_mt1  = 100.*OBS_HM[T=@ngd] / ( OBS_HM[T=@ngd] + OBS_HM[T=@nbd] )
	Let coverage_mt   IF coverage_mt1 GE 100. THEN 99.99 ELSE coverage_mt1
	
	!-- Couverture minimale pour evaluation performances
	!-- ######################################################
	Let/BAD=999.  OBS_HM    IF coverage_mt GE ${COVMIN} THEN OBS_H1   ELSE 999.
	Let/BAD=999.  INNOV_HM  IF coverage_mt GE ${COVMIN} THEN INNOV_H1 ELSE 999.
			
	!!-- Indices performances
	!! ######################
	Let obs_mt     = OBS_HM[ T=@AVE ]
	Let prev       = OBS_HM - INNOV_HM
	Let prev_mt    = prev[ T=@AVE ]
	 
	!! Stats difference OBS-MODELE
	Let innov_mt   = INNOV_HM[ T=@AVE ]
	Let innov_maxt = INNOV_HM[ T=@MAX ]
	Let innov_mint = INNOV_HM[ T=@MIN ]
		
        !! Biais systematique
	Let innov_mima= innov_mint*innov_maxt
	Let recov_tmp1 = IF ( innov_mima NE 999.) THEN 0 ELSE 999.
	Let recov_tmp2 = IF ( innov_mint LT 0. AND innov_mima GT 0. ) THEN (0.75) ELSE recov_tmp1
	Let bias = IF ( innov_mima NE 999.  AND innov_mint GT 0. AND innov_mima GT 0. ) THEN (-0.75) ELSE recov_tmp2
	
	!! Ecart RMS OBS - MODELE
	Let  i2 = INNOV_HM^2.
	Let  RMSe_mt  = i2[ T=@AVE ]^0.5

	!! Ecart RMS centre OBS - MODELE
	Let  ic2 = (INNOV_HM - innov_mt)^2.
	Let  RMSe_cent_mt  = ic2[ T=@AVE ]^0.5

	!! Variance
	Let obs_std_mt     = OBS_HM[ T=@VAR ]^0.5
	Let prev_std_mt    = prev[ T=@VAR ]^0.5
	Let d_std	   = obs_std_mt - prev_std_mt
		
	!! FIGURES
	!! ####################
	
	!! Parametrage
	Define Symbol OUT ${TYPE_FIG} ! PS ou GIF
	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN Set Mode Metafile
	Define Symbol MODEL ${MODEL_TYPE}       
	
	!! Time mean : Data coverage
	!! ##############################
	Shade /Lev=(0,100,10)  coverage_mt
   	ppl Title "${OBS_TYPE}  coverage ($DATE1:$DATE2)"
   	ppl Shade
	Contour /Nolab /Over /Lev=(30,100,30) coverage_mt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG0}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG0}.gif"
	ENDIF	
	
	!! Time mean : SST obs
	!! ##############################
   	!!Shade /Nolab /Lev=(10)(15,30,1)  /Set obs_mt
	Shade /Lev="(8)(12,14,.2) (18)"  obs_mt
   	ppl Title "[ ${OBS_TYPE}, ${MODEL_TYPE} ] : mean SST ${OBS_TYPE} [(^oC, $DATE1:$DATE2)"
   	ppl Shade
	Contour /Nolab /Over /Lev=(15,30,2) obs_mt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG1}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG1}.gif"
	ENDIF

	!! Time mean : SST modele
	!! ##############################
   	!!Shade /Nolab /Lev=(10)(15,30,1) /Set prev_mt
	Shade /Lev="(8)(12,14,.2) (18)"  prev_mt
   	ppl Title "[ ${OBS_TYPE}, ${MODEL_TYPE} ] : mean SST (\$model) [(^oC, $DATE1:$DATE2)"
   	ppl Shade
	Contour /Nolab /Over /Lev=(15,30,2) prev_mt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG2}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG2}.gif"
	ENDIF

	!! Time mean Innovation
	!! ######################################
   	Shade /Lev=(-5)(-3,3,0.5)(5) innov_mt
   	ppl Title "Bias = ${OBS_TYPE} - ${MODEL_TYPE} (^oC, $DATE1:$DATE2)"
   	ppl Shade
   	Contour /Nolab /Over /Lev=(-3)(-2)(-1)(1)(2)(3) innov_mt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG3}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG3}.gif"
	ENDIF

	!! Time mean  RMS
	!! ################################        
	Shade /Lev=(0.,3.5,0.5)(5) RMSe_mt
   	ppl Title "Full RMS difference [ ${OBS_TYPE}, ${MODEL_TYPE} ] (^oC, $DATE1:$DATE2)"
   	ppl Shade
   	Contour /Nolab /Over /Lev=(1)(3)     RMSe_mt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
        !! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG4}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG4}.gif"
	ENDIF

	!! Time mean  centered RMS
	!! ################################        
   	!!Shade /Nolab /Lev=( 0, 3., 0.5 )(5)(10) /Set RMSe_cent_mt
	Shade /Lev=(0.,3.5,0.5)(5) RMSe_cent_mt
   	ppl Title "Centered RMS difference [ ${OBS_TYPE}, ${MODEL_TYPE} ] (^oC, $DATE1:$DATE2)"
   	ppl Shade
   	Contour /Nolab /Over /Lev=(1)(3)     RMSe_cent_mt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
        !! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG4}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG4bis}.gif"
	ENDIF

	!! Innov min
	!! #########       
   	!!Shade /Nolab /Lev=(-10)(-5)( -3, 3, 0.5 )(-5)(10) /Set innov_mint
	Shade /Lev=(-5,5,1) innov_mint
   	ppl Title "Minimum Bias = Min ( ${OBS_TYPE} - ${MODEL_TYPE} ) (^oC, $DATE1:$DATE2)"
   	ppl Shade
   	Contour /Nolab /Over /Lev=(-3)(-1)(1)(3) innov_mint
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG5}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG5}.gif"
	ENDIF

	!! Innov max
	!! #########
   	!!Shade /Nolab /Lev=(-10)(-5)( -3, 3, 0.5 )(-5)(10) /Set innov_maxt
	Shade /Lev=(-5,5,1) innov_maxt
   	ppl Title "Maximum Bias = Max ( ${OBS_TYPE} - ${MODEL_TYPE} ) (^oC, $DATE1:$DATE2)"
   	ppl Shade
   	Contour /Nolab /Over /Lev=(-3)(-1)(1)(3) innov_maxt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG6}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG6}.gif"
	ENDIF
	
	!! Bias 
	!! ####	
	shade /lev=(-1)(-0.5)(0.5)(1) bias
 	ppl Title "[ ${OBS_TYPE}, ${MODEL_TYPE} ] : Systematic Bias ($DATE1:$DATE2)"
   	ppl Shade
   	Contour /Nolab /Over /Lev=(-0.75)(0.75) bias
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG7}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG7}.gif"
	ENDIF	
	
	!! (variance OBS)^0.5 
	!! ####	
	shade /lev=(0,1,0.2)(1.5,3,0.5)(5) obs_std_mt
 	ppl Title "(SST ${OBS_TYPE} Variance  )^{0.5} (^oC, $DATE1:$DATE2)"
   	ppl Shade
   	Contour /Nolab /Over /Lev=(0.5)(1)(1.5)(2) obs_std_mt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG8}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG8}.gif"
	ENDIF
	
	!! (variance PREV)^0.5 
	!! ####	
	shade /lev=(0,1,0.2)(1.5,3,0.5)(5) prev_std_mt
 	ppl Title "(SST ${MODEL_TYPE} Variance)^{0.5} (^oC, $DATE1:$DATE2)"
   	ppl Shade
   	Contour /Nolab /Over /Lev=(0.5)(1)(1.5)(2) prev_std_mt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG9}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG9}.gif"
	ENDIF
	
	!! (variance OBS)^0.5  - (variance PREV)^0.5 
	!! ##########################################
	shade /lev=(-3)(-2)(-1,1,0.25)(2)(3) obs_std_mt - prev_std_mt
 	ppl Title "(SST ${OBS_TYPE} Variance)^{0.5} - (SST ${MODEL_TYPE} Variance)^{0.5} (^oC, $DATE1:$DATE2)"
   	ppl Shade
   	Contour /Nolab /Over /Lev=(-2)(-1)(-0.5)(0.5)(1)(2) obs_std_mt - prev_std_mt
   	Contour /Nolab /Over /Lev=(200)(1000)(4000)  /Color=1 ${BATHY}

	shade/overlay/palette=grey_light.spk/nolabels mask_terre
   	!! Go ${SCRIPT_LAND_DETAIL}

	IF \`STRCMP("(\$OUT)", "PS") EQ 0\` THEN
  		ppl clsplt
  		Sp Fprint -l cps -p portrait -o ${FIG10}.ps metafile.plt
	ENDIF
	IF \`STRCMP("(\$OUT)", "GIF") EQ 0\` THEN
  		Frame /File="${FIG10}.gif"
	ENDIF
	
EOF

ferret -nojnl -gif -script  ${SCRIPT_FERRET}

#-- Nettoyage
#----------------------------------------------------------------------------
rm ${SCRIPT_FERRET}
if [ $TYPE_FIG == "PS" ] ; then
	rm metafile*.plt*
fi	
#
#-- copy dans la directory cible pour html
#\cp $DIR_FIG/*meanT*.gif /home2/sparta/ftp/SOCOM/EVAL_MENOR/.
exit 0
