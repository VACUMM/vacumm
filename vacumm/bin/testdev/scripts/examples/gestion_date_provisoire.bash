
#!/bin/bash
# Fonctions de datation
. libdate.bash
#echo "An: $1"
JJ=$(date2julien $1 $2 $3) # today julian
Y0=${4:0:4}
M0=${4:5:2}
D0=${4:8:2}
#echo "$Y0 $M0 $D0"
J0=$(date2julien $Y0 $M0 $D0)	  # jour julien date debut
echo " --> JJ: $JJ --> J0: $J0"


#setenv / export

