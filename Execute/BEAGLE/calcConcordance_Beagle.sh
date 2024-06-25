#!/bin/bash

if [ $# != 3 ]
then
	echo "Please execute: parallel -j10 -u \"bash calcR2.sh\" $condition $rep < filters "
	exit
fi

condition="$1"
rep="$2"
filter="$3"

wd="$(pwd)"

invs=$(ls $wd/Results_Beagle/$condition/$rep/Inversions)

## For testing
#invs="HsInv1111"
#condition="Standard"
#filter="0.90"

mkdir -p $wd/Results_Beagle/$condition/$rep/Lists/$filter

> $wd/Results_Beagle/$condition/$rep/Lists/$filter/Imp_AFR.conc
> $wd/Results_Beagle/$condition/$rep/Lists/$filter/Imp_EUR.conc
> $wd/Results_Beagle/$condition/$rep/Lists/$filter/Imp_EAS.conc

for inv in $invs
do
	mkdir -p $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter
	chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')
	cat $wd/Results_Beagle/$condition/$rep/Inversions/$inv/ImpGT_AFR $wd/Results_Beagle/$condition/$rep/Inversions/$inv/ImpGT_EUR $wd/Results_Beagle/$condition/$rep/Inversions/$inv/ImpGT_EAS > $wd/Results_Beagle/$condition/$rep/Inversions/$inv/ImpGT 

	Rscript $wd/BEAGLE/makeConcordanceVCF.R $inv $filter $condition $rep
	
	plink --vcf $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/AFR_Exp.vcf --make-bed --out $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/AFR_Exp
	plink --vcf $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/AFR_Imp.vcf --make-bed --out $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/AFR_Imp
	
	plink --bfile $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/AFR_Exp --bmerge $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/AFR_Imp --merge-mode 6 --out $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/AFR
	
	
	plink --vcf $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EUR_Exp.vcf --make-bed --out $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EUR_Exp
	plink --vcf $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EUR_Imp.vcf --make-bed --out $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EUR_Imp
	
	plink --bfile $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EUR_Exp --bmerge $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EUR_Imp --merge-mode 6 --out $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EUR
	
	
	plink --vcf $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EAS_Exp.vcf --make-bed --out $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EAS_Exp
	plink --vcf $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EAS_Imp.vcf --make-bed --out $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EAS_Imp
	
	plink --bfile $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EAS_Exp --bmerge $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EAS_Imp --merge-mode 6 --out $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EAS

#	rm $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/AFR.vcf $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EUR.vcf $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EAS.vcf 

	paste <(echo $inv) <(tail -3 $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/AFR.log | head -n 1 | awk '{print $8}' | sed 's/.$//') >> $wd/Results_Beagle/$condition/$rep/Lists/$filter/Imp_AFR.conc
	paste <(echo $inv) <(tail -3 $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EUR.log | head -n 1 | awk '{print $8}' | sed 's/.$//') >> $wd/Results_Beagle/$condition/$rep/Lists/$filter/Imp_EUR.conc
	paste <(echo $inv) <(tail -3 $wd/Results_Beagle/$condition/$rep/R2/$inv/$filter/EAS.log | head -n 1 | awk '{print $8}' | sed 's/.$//') >> $wd/Results_Beagle/$condition/$rep/Lists/$filter/Imp_EAS.conc

done





