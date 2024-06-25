#!/bin/bash

wd="$(pwd)"

impute_methods="Beagle Minimac4 Impute5"
conditions="Standard Shapeit_Phase"
filters="0.00 0.30 0.50 0.70 0.80 0.90 0.95"
pops="AFR EAS EUR"

for method in $impute_methods
do
	for condition in $conditions
	do
		for filter in $filters
		do
			for pop in $pops
			do
				Rscript $wd/Get_Stats/Get_Indvs_Per_Filt.R $method $pop $condition $filter
			done
		done
	done
done
