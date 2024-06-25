#!/bin/bash

wd="$(pwd)"

condition="$1"
rep="$2"

mkdir -p $wd/Results_Impute2/$condition/$rep/Lists/$filter

parallel -j16 -u "bash $wd/IMPUTE2/calcConcordance_Impute2.sh $condition $rep" < $wd/../VCFs/30X/filters


