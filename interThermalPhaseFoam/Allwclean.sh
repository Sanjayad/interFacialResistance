#!/bin/bash
 
#Clean script for interThermalPhaseChangeFoam
wclean incompressibleTwoPhaseThermalMixture
wclean interThermalPhaseChangeFoam
cd Libraries
source Allwclean.sh
cd ../utilities
wclean interThermalPhaseChangeFoam_Cond
