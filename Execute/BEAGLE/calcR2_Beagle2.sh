#!/bin/bash

## TEST FILE, DO NOT RUN

if [ $# != 2 ]
then
	echo "Please execute: parallel -j16 -u \"bash calcR2.sh\" \$condition < filters "
	exit
fi


condition="$1"
filter="$2"

wd="$(pwd)"

invs=$(ls $wd/Results_Beagle/$condition/Inversions)

## For testing
#invs="HsInv0015 HsInv0030"
#condition="Standard"
#filter="0.90"

mkdir -p $wd/Results_Beagle/$condition/Lists/$filter

> $wd/Results_Beagle/$condition/Lists/$filter/Imp_AFR.r2
> $wd/Results_Beagle/$condition/Lists/$filter/Imp_EUR.r2
> $wd/Results_Beagle/$condition/Lists/$filter/Imp_EAS.r2

for inv in $invs
do
	mkdir -p $wd/Results_Beagle/$condition/R2/$inv/$filter
	chr=$(grep "$inv[[:space:]]" $wd/../VCFs/30X/Common_hg38 | cut -f2 | sed 's/chr//g')
	cat $wd/Results_Beagle/$condition/Inversions/$inv/ImpGT_AFR $wd/Results_Beagle/$condition/Inversions/$inv/ImpGT_EUR $wd/Results_Beagle/$condition/Inversions/$inv/ImpGT_EAS > $wd/Results_Beagle/$condition/Inversions/$inv/ImpGT 

	Rscript $wd/BEAGLE/makeDummyVCF_Beagle.R $inv $filter $condition

	if [ $chr == "X" ]
	then
		plink --vcf $wd/Results_Beagle/$condition/R2/$inv/$filter/AFR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_Beagle/$condition/R2/$inv/$filter/AFR_R2
		plink --vcf $wd/Results_Beagle/$condition/R2/$inv/$filter/EUR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_Beagle/$condition/R2/$inv/$filter/EUR_R2
		plink --vcf $wd/Results_Beagle/$condition/R2/$inv/$filter/EAS.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_Beagle/$condition/R2/$inv/$filter/EAS_R2
	
	else
		plink --vcf $wd/Results_Beagle/$condition/R2/$inv/$filter/AFR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --out $wd/Results_Beagle/$condition/R2/$inv/$filter/AFR_R2
		plink --vcf $wd/Results_Beagle/$condition/R2/$inv/$filter/EUR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --out $wd/Results_Beagle/$condition/R2/$inv/$filter/EUR_R2
		plink --vcf $wd/Results_Beagle/$condition/R2/$inv/$filter/EAS.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_Beagle/$condition/R2/$inv/$filter/EAS_R2
	
	fi

	rm $wd/Results_Beagle/$condition/R2/$inv/$filter/AFR.vcf $wd/Results_Beagle/$condition/R2/$inv/$filter/EUR.vcf $wd/Results_Beagle/$condition/R2/$inv/$filter/EAS.vcf  

	paste <(echo $inv) <(tail -1 $wd/Results_Beagle/$condition/R2/$inv/$filter/AFR_R2.ld | awk '{print $7}') >> $wd/Results_Beagle/$condition/Lists/$filter/Imp_AFR.r2
	paste <(echo $inv) <(tail -1 $wd/Results_Beagle/$condition/R2/$inv/$filter/EUR_R2.ld | awk '{print $7}') >> $wd/Results_Beagle/$condition/Lists/$filter/Imp_EUR.r2
	paste <(echo $inv) <(tail -1 $wd/Results_Beagle/$condition/R2/$inv/$filter/EAS_R2.ld | awk '{print $7}') >> $wd/Results_Beagle/$condition/Lists/$filter/Imp_EAS.r2
		
done





