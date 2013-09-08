#!/bin/bash

######################################
# INIT
######################################
#
SRC=$(pwd)   # La ou installer les sources
PREFIX=$HOME/soft/cdat-vacumm   # La ou installer CDAT et le reste
INSTALL_PYQT4=$(false)
echo EDITEZ MOI && exit
export FC=gfortran
export F90=gfortran
#
export PATH=$PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PREFIX/Externals/lib # hack pour CDAT


######################################
# CDAT  (necessaire)
######################################
cd $SRC
# telechargement svn corrompu var cmor indisponible
#svn export http://www-pcmdi.llnl.gov/svn/repository/cdat/tags/5.2 cdat-5.2
#cd cdat-5.2
# -> on utilise le paquet
wget http://sourceforge.net/projects/cdat/files/Releases/5.2/cdat-5.2.tar.bz2/download
cd cdat-5.2-src
./configure --prefix=$PREFIX #--enable-opendap
make
cd ..

######################################
# SPHINX + extension (necessaire)
######################################
easy_install http://pypi.python.org/packages/2.5/S/Sphinx/Sphinx-1.0-py2.5.egg
easy_install http://pypi.python.org/packages/2.5/s/sphinxcontrib-cheeseshop/sphinxcontrib_cheeseshop-0.2-py2.5.egg

######################################
# CONFIGOBJ (necessaire)
######################################
easy_install http://www.voidspace.org.uk/downloads/configobj-4.6.0.zip

######################################
# PIL (optionel)
######################################
easy_install http://effbot.org/media/downloads/PIL-1.1.7.tar.gz

######################################
# MATPLOTLIB  (necessaire)
######################################
easy_install http://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.0/matplotlib-1.0.0.tar.gz

######################################
# BASEMAP (necessaire)
######################################
BASEMAP_VERSION=1.0
wget http://downloads.sourceforge.net/project/matplotlib/matplotlib-toolkits/basemap-$BASEMAP_VERSION/basemap-$BASEMAP_VERSION.tar.gz
tar xzf basemap-$BASEMAP_VERSION.tar.gz
cd basemap-$BASEMAP_VERSION/geos-3.*/
export GEOS_DIR=$PREFIX/Externals/geos
./configure --prefix=$GEOS_DIR
make install
cd ..
python setup.py install
cd ..

######################################
# PYEXCELERATOR (necessaire)
######################################
PYEXCELERATOR_VERSION=0.6.4.1
wget http://pypi.python.org/packages/source/p/pyExcelerator/pyexcelerator-$PYEXCELERATOR_VERSION.tar.bz2
tar xjf pyexcelerator-$PYEXCELERATOR_VERSION.tar.bz2
cd pyexcelerator-$PYEXCELERATOR_VERSION
python setup.py install
cd ..
# easy_install http://pypi.python.org/packages/source/p/pyExcelerator/pyexcelerator-0.6.4.1.tar.bz2

######################################
# PARAMIKO (necessaire)
######################################
easy_install paramiko


######################################
# SEAWATER (optionel)
######################################
SEAWATER_VERSION=1.0.5
wget http://pypi.python.org/packages/source/s/seawater/seawater-$SEAWATER_VERSION.tar.bz2
tar xjf seawater-$SEAWATER_VERSION.tar.bz2
cd seawater-$SEAWATER_VERSION
python setup.py install
cd ..
# easy_install http://pypi.python.org/packages/source/s/seawater/seawater-1.0.5.tar.bz2


######################################
# PyQt4 (optionel)
######################################
if test $INSTALL_PYQT4 ; then

    # Qt4
    QT4_PACKAGE=qt-everywhere-opensource-src-4.6.3
    wget http://get.qt.nokia.com/qt/source/$QT4_PACKAGE.tar.gz
    tar xzf $QT4_PACKAGE.tar.gz
    cd $QT4_PACKAGE
    ./configure -prefix $PREFIX/Externals/qt
    make && make install
    cd ..
    export QTDIR=$PREFIX/Externals/qt
    export PATH=$QTDIR/bin:$PATH
    export LD_LIBRARY_PATH=$QTDIR/lib:$LD_LIBRARY_PATH

    # SIP
    SIP_PACKAGE=sip-4.10.2
    wget http://www.riverbankcomputing.co.uk/static/Downloads/sip4/$SIP_PACKAGE.tar.gz
    tar xzf $SIP_PACKAGE.tar.gz
    cd $SIP_PACKAGE
    python configure.py
    make && make install
    cd ..

    # PyQt4
    PYQT4_PACKAGE=PyQt-x11-gpl-4.7.3
    wget http://www.riverbankcomputing.co.uk/static/Downloads/PyQt4/$PYQT4_PACKAGE.tar.gz
    tar xzf $PYQT4_PACKAGE.tar.gz
    cd $PYQT4_PACKAGE
    python configure.py
    make && make install 
    cd ..
        

fi

