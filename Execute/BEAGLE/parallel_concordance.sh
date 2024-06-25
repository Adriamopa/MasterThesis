#!/bin/bash

wd="$(pwd)"

condition="$1"
rep="$2"

mkdir -p $wd/Results_Beagle/$condition/$rep/Lists/$filter

parallel -j16 -u "bash $wd/BEAGLE/calcConcordance_Beagle.sh $condition $rep" < $wd/../VCFs/30X/filters


