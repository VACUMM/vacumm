#!/bin/bash
#
# $Id: libdate.bash,v 1.3 2006/06/23 13:39:04 coat Exp $
#
. libTools.bash

function date2julien {
    if [ $# -ne 3 ]; then
	echo "usage: date2julien Y M D"
	return 1
    fi
    YY=$1
    MM=$2
    DD=$3
   
    YY=$( echo $YY | awk '{ printf("%d", $1) }') 
    MM=$( echo $MM | awk '{ printf("%d", $1) }') 
    DD=$( echo $DD | awk '{ printf("%d", $1) }') 

    GGG=1
    if [ $YY -eq 1582 ]; then
	GGG=0
    fi
    if [ $YY -le 1582 ] && [ $MM -lt 10 ]; then
	GGG=0
    fi
    if [ $YY -le 1582 ] && [ $MM -eq 10 ] && [ $DD -lt 5 ]; then 
	GGG=0
    fi
    JD=$((-1 * $((7 * $(($(( $(($MM + 9)) / 12)) + $YY)) / 4))))
    
    S=1
    A=$(($MM - 9))
    if [ $MM -lt 9 ]; then 
	S=-1
	A=$((9 - $MM))
    fi
    J1=$(($YY + $S * $(($A / 7))))
    J1=$((-1 * $(($(( $(($J1 / 100)) + 1)) * 3 / 4))))
    JD=$((JD + $((275 * $MM / 9)) + $DD + $(($GGG * $J1))))
    JD=$((JD + 1721027 + 2 * $GGG + 367 * $YY))
    
    echo $(($JD - 2433283))
    return 0
}
#
#
#
function julien2date {
    if [ $# -ne 1 ]; then
	echo "usage: julien2date JD"
	return 1
    fi
    JD=$1

    Z=$(($JD+2433283))
    
    if [ $Z -lt 2299161 ]; then
	A=$Z
    else
	I=$(echo "(($Z - 1867216.25)/36524.25)" | bc)
	A=$(($Z + 1 + $I - $(($I/4)) ))
    fi
    B=$(($A + 1524))
    C=$(echo "(($B - 122.1)/365.25)" | bc)
    D=$(echo "$C * 365.25 / 1" | bc)
    T=$(echo "(($B - $D)/ 30.6001)" | bc)
    JJ=$(($B - $D - $(echo "(30.6001 * $T) / 1" | bc) ))
    if [ $T -lt 14 ]; then
	MM=$(($T - 1))
    else
	if [ $T == 14 ] || [ $T == 15 ]; then 
	    MM=$(($T - 13 ))
	fi
    fi
    if [ $MM -gt 2 ]; then
	AA=$(($C - 4716 ))
    else
	if [ $MM == 1 ] || [ $MM == 2 ]; then 
	    AA=$(($C - 4715))
	fi
    fi
    MM=$(echo $MM | awk '{ printf("%.02d", $1) }' )
    JJ=$(echo $JJ | awk '{ printf("%.02d", $1) }' )
    echo $AA $MM $JJ
}
