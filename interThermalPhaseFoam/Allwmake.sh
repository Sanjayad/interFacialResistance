#!/bin/bash
 
#Build script for interThermalPhaseChangeFoam
wmake libso incompressibleTwoPhaseThermalMixture
wmake interThermalPhaseChangeFoam
cd Libraries
source Allwmake.sh
cd ../utilities
wmake interThermalPhaseChangeFoam_Cond
#wmake interThermalPhaseChangeFoam_Smoothing
cd ..
