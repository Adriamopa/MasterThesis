#!/bin/bash

wd="$(pwd)"
condition="$1"
vcfs="$2"
rep="$3"

mkdir -p $wd/Results_Minimac4/$condition/$rep
mkdir -p $wd/Results_Minimac4/$condition/$rep/Inversions

parallel -j16 -u "bash $wd/Minimac4/run_minimac.sh $condition $vcfs $rep" < $wd/../VCFs/30X/Invs2Imp




