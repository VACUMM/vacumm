#!/bin/bash
#####################################################################
#
# PURPOSE:
#    Script for batch submission of python scripts.
#    For
#
# USAGE:
#    python.batch myscript.py --myoptions myargs
#
# EXAMPLE:
#    python.batch evalref.py --log-level=INFO ref2
# 
# SETUP:
#   Please edit your environnement setup pour your config below.
#   You maye alternatively set the VACUMM_MODULES environment 
#   variable to load a vacumm environment module 
#   using "module load vacumm".
#   
#
#####################################################################

# Check args
usage(){ echo "   Usage: python.batch myscript.py --myoptions myargs";}
if test "$1" = "-h" -o "$1" = "--help" ; then
    echo "Script for batch submission of python scripts"
    usage
    exit 0
elif ! test -f "$1" ; then
    usage
    exit 1
fi

# Names
jobname=$(echo $1 | cut -d. -f1)
batchfile=.batchfile

# Create the batch file
cat > $batchfile << EOF
#PBS -q sequentiel
#PBS -N $jobname
#PBS -l mem=6gb

echo "Loading python environment"
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!! EDIT HERE YOUR ENVIRONMENT SETUP !!!!!!!!!
if test -n "$VACUMM_MODULES" ; then
    source /usr/share/modules/init/csh
    module use $VACUMM_MODULES
    module load vacumm
fi
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

echo "Changing directory to $(pwd)"
cd $(pwd)

echo "Running python script"
python $1 ${*:2} || exit 1

echo "Done"

EOF

# Submit it
qsub -S /bin/bash $batchfile
