#!/bin/bash

wd="$(pwd)"
rep="$1"


invs=$(ls $wd/Results_ScoreInvHap/Inputs)

## For testing
#invs="HsInv1111"
#condition="Standard"
#filter="0.90"

> $wd/Results_ScoreInvHap/$rep/Imp_AFR.conc
> $wd/Results_ScoreInvHap/$rep/Imp_EUR.conc
> $wd/Results_ScoreInvHap/$rep/Imp_EAS.conc

for inv in $invs
do
	mkdir -p $wd/Results_ScoreInvHap/$rep/R2/$inv
	chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')
	
	cat <(cat $wd/Results_ScoreInvHap/$rep/Outputs/$inv/sIH_output_AFR.csv) <(tail -n+2 $wd/Results_ScoreInvHap/$rep/Outputs/$inv/sIH_output_EUR.csv) <(tail -n+2 $wd/Results_ScoreInvHap/$rep/Outputs/$inv/sIH_output_EAS.csv) > $wd/Results_ScoreInvHap/$rep/Outputs/$inv/sIH_output.csv 

	Rscript $wd/scoreInvHap/makeConcordanceVCF.R $inv $rep
	
	plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_Exp.vcf --make-bed --out $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_Exp
	plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_Imp.vcf --make-bed --out $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_Imp
	
	plink --bfile $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_Exp --bmerge $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_Imp --merge-mode 6 --out $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR
	
	
	plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_Exp.vcf --make-bed --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_Exp
	plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_Imp.vcf --make-bed --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_Imp
	
	plink --bfile $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_Exp --bmerge $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_Imp --merge-mode 6 --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR
	
	
	plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_Exp.vcf --make-bed --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_Exp
	plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_Imp.vcf --make-bed --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_Imp
	
	plink --bfile $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_Exp --bmerge $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_Imp --merge-mode 6 --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS

#	rm $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR.vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR.vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS.vcf 


	paste <(echo $inv) <(if [ $(awk -F'\t' '{print NF; exit}' "$wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_Exp.vcf") -ge 10 ]; then tail -3 $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR.log | head -n 1 | awk '{print $8}' | sed 's/.$//'; else echo "Monomorphic"; fi) >> $wd/Results_ScoreInvHap/$rep/Imp_AFR.conc
	paste <(echo $inv) <(if [ $(awk -F'\t' '{print NF; exit}' "$wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_Exp.vcf") -ge 10 ]; then tail -3 $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR.log | head -n 1 | awk '{print $8}' | sed 's/.$//'; else echo "Monomorphic"; fi) >> $wd/Results_ScoreInvHap/$rep/Imp_EUR.conc
	paste <(echo $inv) <(if [ $(awk -F'\t' '{print NF; exit}' "$wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_Exp.vcf") -ge 10 ]; then tail -3 $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS.log | head -n 1 | awk '{print $8}' | sed 's/.$//'; else echo "Monomorphic"; fi) >> $wd/Results_ScoreInvHap/$rep/Imp_EAS.conc

done





