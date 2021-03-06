#!/bin/bash
source cleanup.sh
rm 0/* DataSummary.dat foam.log
cp A/* 0/
unset FOAM_SIGFPE
blockMesh

#Need to check if makeAxialMesh is available
HAS_MAKEAXIALMESH=`command -v makeAxialMesh`
if [ ! $HAS_MAKEAXIALMESH ]
then
	echo "Need to install makeAxialMesh"
	mkdir -p $HOME/Downloads
	CURDIR=`pwd`
	cd $HOME/Downloads
	svn checkout svn://svn.code.sf.net/p/openfoam-extend/svn/trunk/Breeder_2.0/utilities/mesh/manipulation/MakeAxialMesh
	cd MakeAxialMesh
	./Allwmake
	cd $CURDIR
fi

makeAxialMesh -overwrite
collapseEdges -overwrite
setFields

interThermalPhaseChangeFoam



