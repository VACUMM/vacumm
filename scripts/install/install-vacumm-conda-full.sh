#!/bin/bash

usage() {
    echo "$(basename $0) [-h|--help] [{-p|--prefix} PREFIX] [{-c|--condainstaller} CONDA_INSTALLER]"
}

KERNEL=$(uname -s)
HARDWARE=$(uname -i)
DEFAULT_PREFIX=$HOME/miniconda2
DEFAULT_CONDA_INSTALLER=$HOME/src/Miniconda2-latest-$KERNEL-$HARDWARE.sh

_help() {
    echo "Usage:"
    echo
    usage
    echo
    echo " PREFIX: installation prefix (default: $DEFAULT_PREFIX)"
    echo " CONDA_INSTALLER: conda installer full path (default: $DEFAULT_CONDA_INSTALLER)"
}

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -h|--help)
        _help
        exit 0
        ;;
        -p|--prefix)
        PREFIX="$2"
        shift
        ;;
        -c|--condainstaller)
        CONDA_INSTALLER="$2"
        shift
        ;;
        *)
        echo "Wrong usage"
        usage
        exit 1
        ;;
    esac
    shift
done

# Default values
test -z $PREFIX && PREFIX=$DEFAULT_PREFIX
test -z $CONDA_INSTALLER && CONDA_INSTALLER=$DEFAULT_CONDA_INSTALLER

# Show config
echo "Will install conda with vacumm and dependencies in $PREFIX with installer $CONDA_INSTALLER"
echo "Ok?"
read

# Erase
if test -d "$PREFIX" ; then
    echo "Remove existing installation ?"
    echo "  $PREFIX"
    read
    rm -rf $PREFIX
fi

# Installer
if ! test -f $CONDA_INSTALLER ; then
    echo "Downloading conda installer..."
    mkdir -p $(dirname $CONDA_INSTALLER)
    wget -O $CONDA_INSTALLER https://repo.continuum.io/miniconda/$(basename $CONDA_INSTALLER) || ( echo "Error downloading conda installer" ; exit 1)
fi
chmod +x $CONDA_INSTALLER

# Install conda
$CONDA_INSTALLER -p $PREFIX -b || ( echo "Installation of conda failed" && exit 1 )
export PATH=$PREFIX/bin:$PATH

# Install vacumm and dependencies
conda install -y -c uvcdat uvcdat || ( echo "Installation of uvcdat failed" && exit 1 )
conda install -y -c conda-forge configobj PIL paramiko xlutils seawater pytz cmocean || ( echo "Installation of other dependencies failed" && exit 1 )
conda install -y -c vacumm vacumm || ( echo "Installation of vacumm failed" && exit 1 )

echo
echo "To test it:"
echo "  $PREFIX/bin/python -c 'import vcmq'"

