#!/bin/bash
#!/export/home/services/bin/bash
#
#
# $Id: libTools.bash,v 1.13 2006/09/25 16:01:59 coat Exp $
#
function usage {
    echo "Usage: $*"
    echo
}
#
#
#
function getCachedName {
    FCT=getCachedName
    if [ $# != "1" ]
	then
	usage $FCT filename
	return 1
    fi
    typeset FILE=$1

    echo $(dirname $FILE)/.$(basename $FILE)

    return 0
}
#
#
#
function erreur {
    echo "****ERREUR****"
    echo "$@"
    echo "**************"


    if [ "$HOSTNAME" == "kontan" ]
	then
	echo "$(date -u) $@" | mail -s "*** ERREUR PREVIMER D1 - CAPARMOR ***" previd1@actimar.fr

    else
	echo "$(date -u) $@" | mailx -s "*** ERREUR PREVIMER D1 - CAPARMOR ***" previd1@ifremer.fr
    fi
    return 0
}
#
#
#
function warning {
    echo "****WARNING****"
    echo "$@"
    echo "**************"


    if [ "$HOSTNAME" == "kontan" ]
	then
	echo "$(date -u) $@" | mail -s "*** WARNING PREVIMER D1 - CAPARMOR ***" previd1@actimar.fr

    else
	echo "$(date -u) $@" | mailx -s "*** WARNING PREVIMER D1 - CAPARMOR ***" previd1@ifremer.fr
    fi
    return 0
}
#
#
#
function executer {
    "$@"
    if [ $? -ne 0 ]
	then
	erreur dans "$@"
	exit 2
    fi
    return 0
}
#
#
#
function waitFile {
    FCT="waitFile"
    if [ $# != "2" ]
	then
	usage $FCT file maxtry
	return 1
    fi
    
    typeset FILE=$1
    typeset MAXTRY=$2
    
    NTRY=1
    FLGRUN=1
    while [ "$FLGRUN" -eq "1" ]
      do
      if [ ! -r $FILE ]
	  then
	  if  [ "$NTRY" -gt "$MAXTRY" ]
	      then
	      echo "$(date -u) *** STOP - Too many tries - $FILE not found"
	      return 1
	  else
	      sleep 180
	  fi
	  ((NTRY++))
      else
	  FLGRUN=0
      fi
    done
    
    return 0
}
#
#
#
function waitFileYoungerThan {
    FCT="waitFileYoungerThan"
    if [ $# != "3" ]
	then
	usage $FCT file maxtry maxage
	return 1
    fi
    typeset FILE=$1
    typeset MAXTRY=$2
    typeset MAXAGE=$3

    NTRY=1
    FLGRUN=1
    while [ "$FLGRUN" -eq "1" ]
      do
      if [ ! -r $FILE ]
	  then
	  if  [ "$NTRY" -gt "$MAXTRY" ]
	      then
	      echo "$(date -u) *** STOP - Too many tries - $FILE not found"
	      return 1
	  else
	      sleep 180
	  fi
	  ((NTRY++))
      else
	  AGE=$(fileage $FILE)
	  if [ $AGE -lt $MAXAGE ] && [ $AGE -gt 120 ]
	      then
	      FLGRUN=0
	  elif  [ "$NTRY" -gt "$MAXTRY" ]
	      then
	      echo "$(date -u) *** STOP - Too many tries - $FILE too young"
	      return 1
	  else
	      sleep 180
	  fi
	  ((NTRY++))
      fi
    done

    return 0
}

