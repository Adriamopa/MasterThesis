#!/bin/bash

if [ $# != 3 ]
then
	echo "Please execute: parallel -j16 -u \"bash calcR2.sh\" $condition $rep < filters "
	exit
fi


condition="$1"
rep="$2"
filter="$3"

wd="$(pwd)"

invs=$(ls $wd/Results_Impute2/$condition/$rep/Inversions)

## For testing
#invs="HsInv1111"
#condition="Standard" or condition="Unphased"
#filter="0.90"

mkdir -p $wd/Results_Impute2/$condition/$rep/Lists/$filter

> $wd/Results_Impute2/$condition/$rep/Lists/$filter/Imp_AFR.r2
> $wd/Results_Impute2/$condition/$rep/Lists/$filter/Imp_EUR.r2
> $wd/Results_Impute2/$condition/$rep/Lists/$filter/Imp_EAS.r2

for inv in $invs
do
	mkdir -p $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter
	chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')
	cat $wd/Results_Impute2/$condition/$rep/Inversions/$inv/ImpGT_AFR $wd/Results_Impute2/$condition/$rep/Inversions/$inv/ImpGT_EUR $wd/Results_Impute2/$condition/$rep/Inversions/$inv/ImpGT_EAS > $wd/Results_Impute2/$condition/$rep/Inversions/$inv/ImpGT

	Rscript $wd/IMPUTE2/makeDummyVCF_Impute2.R $inv $filter $condition $rep

	if [ $chr == "X" ]
	then
		plink --vcf $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/AFR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/AFR_R2
		plink --vcf $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EUR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EUR_R2
		plink --vcf $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EAS.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EAS_R2
		
	else
		plink --vcf $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/AFR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --out $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/AFR_R2
		plink --vcf $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EUR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --out $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EUR_R2
		plink --vcf $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EAS.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --out $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EAS_R2

	fi

	rm $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/AFR.vcf $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EUR.vcf $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EAS.vcf 

	paste <(echo $inv) <(tail -1 $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/AFR_R2.ld | awk '{print $7}') >> $wd/Results_Impute2/$condition/$rep/Lists/$filter/Imp_AFR.r2
	paste <(echo $inv) <(tail -1 $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EUR_R2.ld | awk '{print $7}') >> $wd/Results_Impute2/$condition/$rep/Lists/$filter/Imp_EUR.r2
	paste <(echo $inv) <(tail -1 $wd/Results_Impute2/$condition/$rep/R2/$inv/$filter/EAS_R2.ld | awk '{print $7}') >> $wd/Results_Impute2/$condition/$rep/Lists/$filter/Imp_EAS.r2

done





