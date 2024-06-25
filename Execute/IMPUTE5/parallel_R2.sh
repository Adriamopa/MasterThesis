#!/bin/bash

wd="$(pwd)"

condition="$1"
rep="$2"

mkdir -p $wd/Results_Impute5/$condition/$rep/Lists/$filter

parallel -j16 -u "bash $wd/IMPUTE5/calcR2_Impute5.sh $condition $rep" < $wd/../VCFs/30X/filters


