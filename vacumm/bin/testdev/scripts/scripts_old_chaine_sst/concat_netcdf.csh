#!/bin/csh

#################################################################
# Concatenation de fichiers nectdf
################################################################

# get the path for library MKL
source /usr/share/modules/init/csh
module load cmkl/9.1.021
setenv MKL_SERIAL YES
module load netcdf-intel-10/3.6.2

set FILE_IN  = ( $1 )
set FILE_OUT = $2

ncrcat $FILE_IN -o $FILE_OUT

exit 0
