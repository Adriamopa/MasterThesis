#!/bin/bash

wd="$(pwd)"
pops="AFR EUR EAS"
invs=$(cat "$wd/../VCFs/30X/Invs2Imp")

> $wd/Results_Stats/Inversions_SNPs_AFR.txt
> $wd/Results_Stats/Inversions_SNPs_EUR.txt
> $wd/Results_Stats/Inversions_SNPs_EAS.txt

for inv in $invs
do
	for pop in $pops
	do
		num=$(($(bcftools view -H --min-ac 1:minor -M2 -m2 $wd/../VCFs/Inversion_Regions/Unphased/$inv/$pop.vcf.gz | wc -l)-1))
		echo -e "$inv\t$num" >> $wd/Results_Stats/Inversions_SNPs_$pop.txt
	done
done
