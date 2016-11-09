#!/bin/bash
 
#Build script for interThermalPhaseChangeFoam
wmake libso incompressibleTwoPhaseThermalMixture
wmake interThermalPhaseChangeFoam
cd Libraries
source Allwmake.sh
cd ../
wmake utilities/interTherMalPhaseChangeFoam_Cond/
