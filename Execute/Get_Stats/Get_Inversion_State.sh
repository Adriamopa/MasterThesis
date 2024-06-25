#!/bin/bash

wd="$(pwd)"
pops="AFR EUR EAS"
invs=$(cat "$wd/../VCFs/30X/Invs2Imp")

> $wd/Results_Stats/Inversions_State_AFR.txt
> $wd/Results_Stats/Inversions_State_EUR.txt
> $wd/Results_Stats/Inversions_State_EAS.txt

for inv in $invs
do
	for pop in $pops
	do
		state=$(bcftools view -H $wd/../VCFs/Inversion_Regions/Unphased/$inv/$pop.vcf.gz | grep $inv | cut -f 10- | sed 's_/_|_g' | tr "\t" "\n" | sort | uniq -c | wc -l)
		if [ $state -gt 1 ]
		then
			echo -e "$inv\tPolymorphic" >> $wd/Results_Stats/Inversions_State_$pop.txt
		elif [ $state -eq 1 ]
		then
			echo -e "$inv\tMonomorphic" >> $wd/Results_Stats/Inversions_State_$pop.txt
		else
			echo -e "Error in Inversion $inv and Population $Pop!!\n"
		fi
	done
done
