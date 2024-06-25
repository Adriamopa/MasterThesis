#!/bin/bash

wd="$(pwd)"

#inv="$1"


## Test
inv="HsInv0015" 

mkdir -p $wd/../VCFs/Inversion_Regions/Impute2_Phase
mkdir -p $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/

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
cat <(cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -vE "ND$|NA$|Del/|/Del" | cut -f1) <(cut -f1 $wd/../VCFs/30X/Panel30x) | sort | uniq -d > $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Samples
cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -f $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Samples | sort -k1 > $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/InvGenotypes.txt
samps=$(cat $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Samples)

## The VCF
bcftools view -r chr$chr:$BP1-$BP2 -s $(echo $samps | tr " " ",") --min-ac 2:minor -M2 -m2 -Ov $vcf | bcftools norm -d all > $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Region.vcf

bcftools annotate --rename-chrs $wd/../VCFs/30X/chr_names.txt  $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Region.vcf -Ov -o $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/RegionMod.vcf	

std_line=$(echo -e "$chr\t$Midpoint\t$inv\tA\tG\t100\tPASS\t.\tGT")
genos=$(cut -f2 $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/InvGenotypes.txt)

paste <(echo $std_line | sed 's/ /\t/g') <(echo $genos | sed 's/ /\t/g' | sed 's/Std/0/g; s/Inv/1/g') >> $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/RegionMod.vcf
vcf-sort $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/RegionMod.vcf > $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/RegionSorted.vcf

rm $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Samples $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/InvGenotypes.txt $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Region.vcf $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/RegionMod.vcf

$wd/../Programes/vcf2impute_gen.pl -vcf $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/RegionSorted.vcf -gen $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/RegionSorted.gen


if [ $chr == "X" ]
then
	Rscript $wd/Extract_regions_pops/MakeForX_Impute2.R $inv
	impute2 \
	-prephase_g \
	-chrX \
	-buffer 0 kb \
	-m $wd/../VCFs/30X/Maps/chrX_hg38_genmap.txt \
	-g $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/RegionSorted.gen.gz \
	-sample_g $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Reference.sample \
	-int $BP1 $BP2 \
	-Ne 20000 \
	-o $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Out
else
	impute2 \
	-prephase_g \
	-g $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/RegionSorted.gen.gz \
	-o $wd/../VCFs/Inversion_Regions/Impute2_Phase/$inv/Out \
	-buffer 0 kb \
	-m $wd/../VCFs/30X/Maps/chrX_hg38_genmap.txt \
	-int $BP1 $BP2 \
	-Ne 20000
fi			
