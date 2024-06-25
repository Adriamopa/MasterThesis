#!/bin/bash

wd="$(pwd)"

condition="$1" #"Standard"
vcfs="$2" #"Standard"
rep="$3"

mkdir -p $wd/Results_Impute5/$condition/$rep/
mkdir -p $wd/Results_Impute5/$condition/$rep/Inversions

parallel -j16 -u "bash $wd/IMPUTE5/run_combined_Impute5.sh $condition $vcfs $rep" < $wd/../VCFs/30X/Invs2Imp
