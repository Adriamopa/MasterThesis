#!/bin/bash

wd="$(pwd)"

### Extract the VCF from adding our genotypes data

inv="$1"

## For testing
#inv="HsInv0015"


pops="AFR EUR EAS"

mkdir -p $wd/Results_ScoreInvHap/
mkdir -p $wd/Results_ScoreInvHap/Inputs/$inv

for pop in $pops
do

	chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')
	BP1=$(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f5)-200000))
	BP2=$(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f8)+200000))
	Midpoint=$(($(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f6)+$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f7)))/2))

	if [ $chr == "X" ]
	then
		vcf="$wd/../VCFs/30X/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
	else
		vcf="$wd/../VCFs/30X/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
	fi


	nline=$(head -1 $wd/../VCFs/30X/GTypesINVs.csv | sed 's/\t/\n/g' | cat -n | grep "$inv$" | awk '{print $1}')

	## For the inclusion of genotypes
	cat <(cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -vE "ND$|NA$|Del/|/Del" | cut -f1) \
		<(grep $pop $wd/../VCFs/30X/Panel30x | cut -f1) | sort | uniq -d > $wd/Results_ScoreInvHap/Inputs/$inv/Samples
	cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -f $wd/Results_ScoreInvHap/Inputs/$inv/Samples | sort -k1 > $wd/Results_ScoreInvHap/Inputs/$inv/InvGenotypes.txt
	samps=$(cat $wd/Results_ScoreInvHap/Inputs/$inv/Samples)

	## The VCF
	bcftools view -r chr$chr:$BP1-$BP2 -s $(echo $samps | tr " " ",") --min-ac 2:minor -M2 -m2 -Ov $vcf | bcftools norm -d all > $wd/Results_ScoreInvHap/Inputs/$inv/Region.vcf

	bcftools annotate --rename-chrs $wd/../VCFs/30X/chr_names.txt  $wd/Results_ScoreInvHap/Inputs/$inv/Region.vcf -Ov -o $wd/Results_ScoreInvHap/Inputs/$inv/RegionMod.vcf

	## Needed convertion of missing variant names to a CHROM_POSITION_REF_ALT notation
	## Also for diploidification of chrX males
	grep "##" $wd/Results_ScoreInvHap/Inputs/$inv/RegionMod.vcf > $wd/Results_ScoreInvHap/Inputs/$inv/Header
	grep -v "##" $wd/Results_ScoreInvHap/Inputs/$inv/RegionMod.vcf > $wd/Results_ScoreInvHap/Inputs/$inv/Body1

	Rscript $wd/scoreInvHap/CompleteIDs_ScoreInvHap.R $wd/Results_ScoreInvHap/Inputs/$inv $pop

	if [ $chr == "X" ]
	then
		cat $wd/Results_ScoreInvHap/Inputs/$inv/Header $wd/Results_ScoreInvHap/Inputs/$inv/Body2_dip > $wd/Results_ScoreInvHap/Inputs/$inv/RegionMod.vcf
	else
		cat $wd/Results_ScoreInvHap/Inputs/$inv/Header $wd/Results_ScoreInvHap/Inputs/$inv/Body2 > $wd/Results_ScoreInvHap/Inputs/$inv/RegionMod.vcf
	fi
	

	## Bed done as input
	plink --vcf $wd/Results_ScoreInvHap/Inputs/$inv/RegionMod.vcf --make-bed --out $wd/Results_ScoreInvHap/Inputs/$inv/Region_$pop


	## Add genotypes to the normal (males as haploid for X) VCF

	#cat $wd/Inputs/$inv/Header $wd/Inputs/$inv/Body2 > $wd/Inputs/$inv/Region.vcf
	std_line=$(echo -e "$chr\t$Midpoint\t$inv\tStd\tInv\t100\tPASS\t.\tGT")
	genos=$(cut -f2 $wd/Results_ScoreInvHap/Inputs/$inv/InvGenotypes.txt)

	paste <(echo $std_line | sed 's/ /\t/g') <(echo $genos | sed 's/ /\t/g' | sed 's/Std/0/g; s/Inv/1/g') >> $wd/Results_ScoreInvHap/Inputs/$inv/RegionMod.vcf
	vcf-sort $wd/Results_ScoreInvHap/Inputs/$inv/RegionMod.vcf > $wd/Results_ScoreInvHap/Inputs/$inv/RegionSorted.vcf

	plink -vcf $wd/Results_ScoreInvHap/Inputs/$inv/RegionSorted.vcf \
		--r2 --ld-snp $inv \
		--ld-window 999999 \
		--ld-window-kb 9999999 \
		--ld-window-r2 0 \
		--out $wd/Results_ScoreInvHap/Inputs/$inv/Region

	## R2 input
	awk '{print $6"\t"$7}' $wd/Results_ScoreInvHap/Inputs/$inv/Region.ld > $wd/Results_ScoreInvHap/Inputs/$inv/R2_$pop.input

	grep -v "##" $wd/Results_ScoreInvHap/Inputs/$inv/RegionSorted.vcf > $wd/Results_ScoreInvHap/Inputs/$inv/Region_$pop.csv

	rm $wd/Results_ScoreInvHap/Inputs/$inv/Samples \
		$wd/Results_ScoreInvHap/Inputs/$inv/InvGenotypes.txt \
		$wd/Results_ScoreInvHap/Inputs/$inv/Region.vcf \
		$wd/Results_ScoreInvHap/Inputs/$inv/RegionMod.vcf \
		$wd/Results_ScoreInvHap/Inputs/$inv/RegionSorted.vcf \
		$wd/Results_ScoreInvHap/Inputs/$inv/Header \
		$wd/Results_ScoreInvHap/Inputs/$inv/Body1 \
		$wd/Results_ScoreInvHap/Inputs/$inv/Body2 \
		$wd/Results_ScoreInvHap/Inputs/$inv/Body2_dip

done









