#!/bin/bash
################################################################
#
# PURPOSE:
#   
#   Update svn doc branch and rebuild doc if needed.
#   This script is useful for crontab execution.
#
#   WARNING: 
#       If does not check local modifications
#       and must never do it.
#   
# SETUP:
#   
#   Set the VACUMM_SVN environment variable to define
#   where your local copy of the svn tree points to.
#   Arguments are passed to "make"
#
# EXAMPLE OF CRONTAB JOB:
#   
#   File update_vacumm_doc.sh:
#
#       setenv VACUMM_SVN "$HOME/VACUMM/svn"
#       $VACUMM_SVN/bin/update_doc.sh $*
#
# EXAMPLE OF CRONTAB LINE:
#   
#   05,35 * * * * /bin/bash $HOME/bin/update_vacumm_doc.sh html
#
################################################################

# Variable not defined
test -z $VACUMM_SVN && exit

# Directory not found
test -d $VACUMM_SVN || exit


# Svn update
cd $VACUMM_SVN/doc/sphinx
nlines=$(svn up | wc -l)
test $(svn up | wc -l) -gt 1 && make $*
