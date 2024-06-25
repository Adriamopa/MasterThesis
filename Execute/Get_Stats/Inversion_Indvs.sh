#!/bin/bash

wd="$(pwd)"

echo -e "INV\tAFR.Num.Samples\tEAS.Num.Samples\tEUR.Num.Samples" > $wd/Results_Stats/Inv_Inds.txt
while IFS= read -r line; do
        samples_AFR=$(zcat $wd/../VCFs/Inversion_Regions/Standard/$line/AFR.vcf.gz | grep "#CHROM" | cut -f10-)
        samples_EAS=$(zcat $wd/../VCFs/Inversion_Regions/Standard/$line/EAS.vcf.gz | grep "#CHROM" | cut -f10-)
        samples_EUR=$(zcat $wd/../VCFs/Inversion_Regions/Standard/$line/EUR.vcf.gz | grep "#CHROM" | cut -f10-)
        num_samples_AFR=$(echo $samples_AFR | tr " " "\n" | wc -l)
        num_samples_EAS=$(echo $samples_EAS | tr " " "\n" | wc -l)
        num_samples_EUR=$(echo $samples_EUR | tr " " "\n" | wc -l)
        echo -e "$line\t$num_samples_AFR\t$num_samples_EAS\t$num_samples_EUR" >> $wd/Results_Stats/Inv_Inds.txt
done < $wd/../VCFs/30X/Invs2Imp

