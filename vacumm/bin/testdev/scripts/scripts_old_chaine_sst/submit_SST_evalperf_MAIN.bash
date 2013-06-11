#!/bin/bash
#**********************************************************************************
#* MOD : submit_SST_MAIN.bash
#* OBJ : Evaluation des performances modele en SST 
#* EXE : submit_evalperf.bash sur CAPARMOR
#* CRE : P. Craneguy (Actimar) - Projet PRECOC (ANR)
#* VER : 1.0 (25/01/2007)
#**********************************************************************************

# Fonctions de datation
. libdate.bash

#--- PARAMETRAGE DES ACTIONS
#--------------------------------------------------------------------------------
#-- Actions et Environnement de SST
. submit_SST_evalperf_par.bash # parametres des actions
. submit_SST_evalperf_env.bash # environnement de SST

cd ${HDIR}
if   [ ! -d ${DIR_SST} ] ; then  
	mkdir ${DIR_SST} 
fi
if   [ ! -d ${DIR_FIG} ] ; then
	mkdir ${DIR_FIG}
fi

echo "%%--------------------------------------------------------------------------"
echo "%%--------------------------------------------------------------------------"
echo "%%-  FORMATION DES FICHIERS DE SST SAT (NAR & L3P supercollated) 	       -"
echo "%%-          POUR EVALUATION DES PERF DE MARS (MENOR)		       -"
echo "%%--------------------------------------------------------------------------"
echo "%%-	scripts FORM_SST_sat_main.bash & submit_SST_prev_FORM_MAIN.bash   -"		
echo "%%--------------------------------------------------------------------------"

	# -- SST satellite (NAR & L3P)
FILESORT_NAR=$DIR_SST/${FOUT_NAR}
FILESORT_L3P=$DIR_SST/${FOUT_L3P}
DIR_SST_NAR=${DIR_SST}_${SST_A} # Repertoire de stockage des fichiers non-concatenes de SST NAR
DIR_SST_L3P=${DIR_SST}_${SST_B} # Repertoire de stockage des fichiers non-concatenes de SST L3P
##
##mkdir ${DIR_SST_NAR} ${DIR_SST_L3P}
##
if   [ ! -d ${DIR_SST_NAR} ] ; then
	mkdir ${DIR_SST_NAR}
fi

if [ "$lance_form_data_sat" == "OUI" ] ; then
   FORM_SST_sat_main.bash  $HDIR ${DIR_SST_NAR} ${DIR_SST_L3P} $FILESORT_NAR $FILESORT_L3P $DATE_FIN $nSSTdayfile
else
   echo " Pas de formation des SST NAR & L3P"
fi

	# -- SST prevision modele (MENOR)
FILESORT_PREV_MARS=$DIR_SST/${FOUT_PREV_MARS}

	# -- SST prevision modele (MFS)
#FILESORT_PREV_MFS=$DIR_SST/${FOUT_PREV_MFS}
# modif PG
#FILESORT_PREV_MFS=$DIR_SST/${FOUT_PREV_MFS}

if [ "$lance_form_data_prev" == "OUI" ] ; then
   FORM_SST_prev.bash ${HDIR} ${DATE_FIN} ${nSSTdayfile} ${FILESORT_PREV_MARS} 
   ##
   ##FORM_SST_mfs.bash  ${HDIR} ${DATE_FIN} ${nSSTdayfile} ${MFS_TYPE} ${FILESORT_PREV_MFS} 
   # modif PG
   #FORM_SST_mfs.bash  ${HDIR} ${DATE_FIN} ${nSSTdayfile} ${MFS_TYPE} ${FILESORT_PREV_MFS}
   ##
else
   echo " Pas de formation des SST prevues MENOR & MFS"
fi

echo "%%--------------------------------------------------------------------------"
echo "%%--------------------------------------------------------------------------"
echo "%%-  		EVALUATION DES PERFORMANCES EN SST   	       	       -"
echo "%%--------------------------------------------------------------------------"
echo "%%-	scripts submit_SST_NARmodel_collocated.sh & submit_SST_evalperf_fig.sh  -"		
echo "%%--------------------------------------------------------------------------"


#if [ -f $HDIR/$FILESORT_PREV_MARS ] && [ -f $HDIR/$FILESORT_PREV_MFS ] && [ -f $HDIR/$FILESORT_NAR ]
if [ -f $HDIR/$FILESORT_PREV_MARS ] && [ -f $HDIR/$FILESORT_NAR ]
 then

	echo "#--- PERFORMANCE PAR RAPPORT A SST NAR"
	echo "#-------------------------------------------------------------------------"

	#-- parametrisation des scripts
	#	-- modele MARS
	FILE_PREV=$FILESORT_PREV_MARS
	VAR_PREV=$VAR_SST_PREV
	MVAL_PREV=$MVAL_SST_PREV
	VMIN_PREV=$VMIN_SST_PREV
	VAR_BATHY=$VAR_HO_PREV
	MVAL_BATHY=$MVAL_H0_PREV
	PAR_PREV="$FILE_PREV $VAR_PREV $MVAL_PREV $VMIN_PREV $VAR_BATHY $MVAL_BATHY"

	#	-- modele MFS
        FILE_PREV_MFS=$FILESORT_PREV_MFS
	VAR_PREV_MFS=$MFS_VAR_TEMP
	PAR_PREV_MFS="$FILE_PREV_MFS $VAR_PREV_MFS $MVAL_PREV $VMIN_PREV $VAR_BATHY $MVAL_BATHY"

	#	-- OBS
	FILE_OBS=$FILESORT_NAR
	VAR_OBS=$VAR_SST_NAR
	UNIT_OBS=$UNIT_SST_NAR
	VAR_LEV=$VAR_LEV_SST_NAR
	LEV_MIN=$LEV_MIN_COLnarprev
	PAR_OBS="$FILE_OBS $VAR_OBS $UNIT_OBS $VAR_LEV $LEV_MIN"

	#	-- OUT
	#		-- modele MARS
	FILEOUT=${FILEOUT_COLnarprev_NAR}
	VAR_OUT_PREV=${VAR_OUT_COLnarprev_prev}
	VAR_OUT_OBS=${VAR_OUT_COLnarprev_OBS}
	VAR_OUT_BATHY=${VAR_OUT_COLnarprev_bathy}
	MVAL_OUT=$MVAL_OUT
	PAR_OUT="$FILEOUT $VAR_OUT_PREV $VAR_OUT_OBS $VAR_OUT_BATHY $MVAL_OUT"
	
	#		-- modele MFS
        FILEOUT_MFS=${FILEOUT_COLnarmfs_NAR}
	PAR_OUT_MFS="$FILEOUT_MFS $VAR_OUT_PREV $VAR_OUT_OBS $VAR_OUT_BATHY $MVAL_OUT"
	
	#		-- parametres
	DATE_END=${DATE_FIN}
	DUREE=${nSSTdayfile}
	HEURE=$HOUR
	D_HEURE=${D_HOUR}
	PAR_TIME="$DATE_END $DUREE $HEURE $D_HEURE"

	OBS_TYPE=$SST_A
	MODEL_TYPE=MARS
	SCRIPT_LAND_DETAIL=${SCRIPT_LAND_DETAIL}
	PAR_1="$FILEOUT   $VAR_OUT_PREV $VAR_OUT_OBS $VAR_OUT_BATHY"
	PAR_2="$OBS_TYPE $MODEL_TYPE $DATE_END $DUREE $SCRIPT_LAND_DETAIL"
	PAR_3="$DIR_FIG $FIG_FORMAT"
	
	if [ "$lance_eval_perf" == "OUI" ]
	 then
	   # Bathy MARS projettée sur grille NAR
	   MODEL_TYPE=MARS
	   sh COLLOC_NAR_model.sh $PAR_PREV $PAR_OBS $PAR_OUT $PAR_TIME $MODEL_TYPE

	   # Bathy MFS n est pas traitee lors de la collocalisation
	   ##MODEL_TYPE=MFS
	   ##sh COLLOC_NAR_model.sh $PAR_PREV_MFS $PAR_OBS $PAR_OUT_MFS $PAR_TIME $MODEL_TYPE

	   n_zone_m1=$(( $n_zone - 1 ))
	   n=0
	   while [ $n -le $n_zone_m1 ]
	   do	
		echo "--------------- ZONE :" ${ZONE[$n]}
		HLIM=$HMAX_TOT
		if [ "${ZONE[$n]}" == "GOL" ]
		then
			HLIM=$HMAX_SHELF
			
		else
			HLIM=$HMAX_TOT
		fi
	  	#-- moyenne temporelle par zone
	   	sh SST_evalperf_fig_meanT.sh  $PAR_1 $PAR_2 $PAR_3  ${ZONE[$n]}  $HLIM $COV_SST_MIN_MEANT
	   	#-- evolution temporelle par zone
	   	sh SST_evalperf_fig_evolT.sh  $PAR_1 $PAR_2 $PAR_3  ${ZONE[$n]}  $HLIM $COV_SST_MIN_EVOLT
	   	n=$(( $n + 1 ))
	   done	   
	   # rm $FILEOUT
	else 
	    echo " Pas d evaluation des performances model avec SST NAR"
	fi
else
	#echo "--> Le(s) fichier(s) de SST NAR ou (et) Prevision (MARS ou(et) MFS) n existe(nt) pas"
	echo "--> Le(s) fichier(s) de SST NAR ou (et) Prevision MARS n existe(nt) pas"
	echo " SST Mars : " $HDIR/$FILESORT_PREV_MARS 
	echo " SST NAR  : " $HDIR/$FILESORT_NAR
fi
exit 0





