#!/bin/bash

wd="$(pwd)"

condition="$1"
vcfs="$2"
rep="$3"

mkdir -p $wd/Results_Impute2
mkdir -p $wd/Results_Impute2/Inversions/$condition

parallel -j16 -u "bash $wd/run_Impute2.sh $condition $vcfs $rep" < $wd/../VCFs/30X/Invs2Imp
