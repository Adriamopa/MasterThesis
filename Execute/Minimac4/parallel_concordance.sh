#!/bin/bash

wd="$(pwd)"

condition="$1"
rep="$2"

mkdir -p $wd/Results_Minimac4/$condition/$rep/Lists/$filter

parallel -j16 -u "bash $wd/Minimac4/calcConcordance_Minimac4.sh $condition $rep" < $wd/../VCFs/30X/filters


