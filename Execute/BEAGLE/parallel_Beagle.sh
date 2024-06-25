#!/bin/bash

wd="$(pwd)"
condition="$1"
vcfs="$2"
rep="$3"

mkdir -p $wd/Results_Beagle/$condition/$rep
mkdir -p $wd/Results_Beagle/$condition/$rep/Inversions

parallel -j16 -u "bash $wd/BEAGLE/run_Beagle.sh $condition $vcfs $rep" < $wd/../VCFs/30X/Invs2Imp




