#!/bin/bash

wd=$(pwd)
invs=$(ls Results_ScoreInvHap/Inputs)
rep="$1"
for inv in $invs
do
	pops="AFR EUR EAS"
	mkdir -p $wd/Results_ScoreInvHap/$rep/Outputs/$inv
	echo "$inv"
	for pop in $pops
	do
		Rscript $wd/scoreInvHap/sIH.R $inv $pop $rep
	done
done



