#!/bin/bash

wd=$(pwd)

condition="$1"
rep="$2"
filters="0.50 0.70 0.80 0.90 0.95"
mkdir -p $wd/Results_Impute2/$condition/$rep/Lists
#parallel -j16 -u "bash $wd/IMPUTE2/calcR2_Impute2.sh $condition $rep" < $wd/../VCFs/30X/filters

for filter in $filters
do
	bash $wd/IMPUTE2/calcR2_Impute2.sh $condition $rep $filter
done
