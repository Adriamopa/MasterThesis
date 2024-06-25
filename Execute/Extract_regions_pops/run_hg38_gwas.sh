#!/bin/bash

wd="$(pwd)"

inv="$1"


## Test
# condition="Standard" or condition="Unphased"
# inv="HsInv0015" 

mkdir -p $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv
mkdir -p $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/

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
cat <(cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -vE "ND$|NA$|Del/|/Del" | cut -f1) <(cut -f1 $wd/../VCFs/30X/Panel30x) | sort | uniq -d > $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Samples
cut -f1,$nline $wd/../VCFs/30X/GTypesINVs.csv | grep -f $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Samples | sort -k1 > $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/InvGenotypes.txt
samps=$(cat $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Samples)

## The VCF
bcftools view -r chr$chr:$BP1-$BP2 -s $(echo $samps | tr " " ",") --min-ac 2:minor -M2 -m2 -Ov $vcf | bcftools norm -d all > $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Region.vcf

bcftools annotate --rename-chrs $wd/../VCFs/30X/chr_names.txt  $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Region.vcf -Ov -o $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/RegionMod.vcf	

std_line=$(echo -e "$chr\t$Midpoint\t$inv\tA\tG\t100\tPASS\t.\tGT")
genos=$(cut -f2 $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/InvGenotypes.txt)

paste <(echo $std_line | sed 's/ /\t/g') <(echo $genos | sed 's/ /\t/g' | sed 's/Std/0/g; s/Inv/1/g') >> $wd/../VCFs/GWAS_VCFs/Input_VCFs/$/$inv/RegionMod.vcf
vcf-sort $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/RegionMod.vcf > $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/RegionSorted.vcf

#cat $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/RegionSorted.vcf | sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' > $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/RegionFilled.vcf

# rm $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Samples $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/InvGenotypes.txt $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Region.vcf $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/RegionMod.vcf $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/RegionSorted.vcf

java -jar $wd/../Programes/imputation_methods/beagle.22Jul22.46e.jar \
	gt=$wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/RegionSorted.vcf \
	out=$wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Output \
	map=$wd/../VCFs/30X/Maps/plink.chr${chr}.GRCh38.map
	
gunzip $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Output.vcf.gz

sed -i -E 's/END=[0-9]+/./g' $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Output.vcf
sed -i -E "/##source/ a ##contig=<ID=${chr},length=10000000000000000>" $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Output.vcf

bgzip $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Output.vcf

tabix -f -p vcf $wd/../VCFs/GWAS_VCFs/Input_VCFs/$inv/Output.vcf.gz
