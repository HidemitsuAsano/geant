#!/bin/sh
echo "Making knucl4 START"
INSTALL_DIR=build
G4INSTALL_DIR=/group/had/knucl/knucl_soft/geant4.10.01.p01_new_install

ANALYZER_DIR=../../Run68/k18ana

if [ ! -e $INSTALL_DIR ]; then
    echo "Make $INSTALL_DIR Directory"
    mkdir $INSTALL_DIR
fi
cd $INSTALL_DIR
source $G4INSTALL_DIR/bin/geant4.sh
cmake -DGeant4_DIR=$G4INSTALL_DIR/lib64/Geant4-10.1.1 ../
echo "Make knucl in $INSTALL_DIR"
make
make -f Makefile.lib
if [ ! -e $ANALYZER_DIR ]; then
    echo "not find $ANALYZER_DIR"
else
    echo "analyzer dir=$ANALYZER_DIR"
    if [ ! -e conf ]; then
	echo "make link to conf directory"
	ln -s $ANALYZER_DIR/conf
    fi
    echo "analyzer dir=$ANALYZER_DIR"
    if [ ! -e param ]; then
	echo "make link to param directory"
	ln -s $ANALYZER_DIR/param
    fi
fi
if [ ! -e MY.list ]; then
    echo "not find Default CS.list ..."
    echo "cp CS.list -> MY.list"
    cp ../CS.list MY.list
fi
if [ -e $ANALYZER_DIR/src ]; then
    if [ ! -e $ANALYZER_DIR/src/lib/KnuclRootData_cc.so ]; then
	echo "make link $ANALYZER_DIR/src/KnuclRootData_cc.so -> KnuclRootData_cc.so"
	cd $ANALYZER_DIR/src/lib
	if [ -e KnuclRootData_cc.so ]; then
	    rm KnuclRootData_cc.so
	fi
	ln -s ../../../geant/knucl4/build/KnuclRootData_cc.so .
	cd ../
	if [ -e KnuclRootData.h ]; then
	    rm KnuclRootData.h
	fi
	ln -s ../../geant/knucl4/build/KnuclRootData.h .
	if [ -e ComCrossSectionTable.hh ]; then
	    rm ComCrossSectionTable.hh
	fi
	ln -s ../../geant/knucl4/build/ComCrossSectionTable.hh .
    fi
fi