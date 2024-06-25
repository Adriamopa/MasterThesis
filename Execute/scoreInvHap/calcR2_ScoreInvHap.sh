#!/bin/bash

wd="$(pwd)"
rep="$1"


> $wd/Results_ScoreInvHap/$rep/Imp_AFR.r2
> $wd/Results_ScoreInvHap/$rep/Imp_EUR.r2
> $wd/Results_ScoreInvHap/$rep/Imp_EAS.r2

pops="AFR EUR EAS"

for inv in $(ls $wd/Results_ScoreInvHap/Inputs)
do
	mkdir -p $wd/Results_ScoreInvHap/$rep/R2/$inv

	cat <(cat $wd/Results_ScoreInvHap/$rep/Outputs/$inv/sIH_output_AFR.csv) <(tail -n+2 $wd/Results_ScoreInvHap/$rep/Outputs/$inv/sIH_output_EUR.csv) <(tail -n+2 $wd/Results_ScoreInvHap/$rep/Outputs/$inv/sIH_output_EAS.csv) > $wd/Results_ScoreInvHap/$rep/Outputs/$inv/sIH_output.csv

	Rscript $wd/scoreInvHap/makeDummyVCF_ScoreInvHap.R $inv $rep

	chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')

	if [ $chr == "X" ]
	then
		plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_R2
		plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_R2
		plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --ld-xchr 1 --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_R2
	else
		plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --out $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_R2
		plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_R2
		plink --vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS.vcf --r2 --ld-window-kb 1000000 --ld-window 999999 --ld-window-r2 0 --out $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_R2
	fi

	rm $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR.vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR.vcf $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS.vcf 

	paste <(echo $inv) <(if [ -f "$wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_R2.ld" ]; then tail -1 $wd/Results_ScoreInvHap/$rep/R2/$inv/AFR_R2.ld | awk '{print $7}'; else echo "R2"; fi) >> $wd/Results_ScoreInvHap/$rep/Imp_AFR.r2
	paste <(echo $inv) <(if [ -f "$wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_R2.ld" ]; then tail -1 $wd/Results_ScoreInvHap/$rep/R2/$inv/EUR_R2.ld | awk '{print $7}'; else echo "R2"; fi) >> $wd/Results_ScoreInvHap/$rep/Imp_EUR.r2
	paste <(echo $inv) <(if [ -f "$wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_R2.ld" ]; then tail -1 $wd/Results_ScoreInvHap/$rep/R2/$inv/EAS_R2.ld | awk '{print $7}'; else echo "R2"; fi) >> $wd/Results_ScoreInvHap/$rep/Imp_EAS.r2

done


