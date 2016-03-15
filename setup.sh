#!/bin/bash

####################### Prefixes #######################

if ! [ -x "$(command -v applgrid-config)" ]; then
    echo "APPLgrid config file not found in path!"
    echo "aborting installation"
    exit -1
fi

if ! [ -x "$(command -v apfel-config)" ]; then
    echo "APFEL config file not found in path!"
    echo "aborting installation"
    exit -1
fi

APPLVER=applgrid-$(applgrid-config --version)
TARGET=$(applgrid-config --incdir)/appl_grid/
PREFIX=$(apfel-config --prefix)

####################### Header installation #######################

echo "Detected APPLgrid version: " $APPLVER" .. supplementing with full headers ... "
wget http://www.hepforge.org/archive/applgrid/$APPLVER.tgz

if [ ! -f "./"$APPLVER".tgz" ]; then
    echo "APPLgrid tgz failed to download!"
    exit -1
fi

tar -xzf ./$APPLVER.tgz
cp "./"$APPLVER"/src/"*.h $TARGET
rm $APPLVER.tgz
rm -rf ./$APPLVER

 
####################### Installation #######################
echo
echo "Headers supplemented, proceeding to installation in APFEL directory: "$PREFIX
echo
autoreconf -i
./configure --prefi=$PREFIX
make && make check && make install

exit 1