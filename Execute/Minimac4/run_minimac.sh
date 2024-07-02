#!/bin/bash

if [ $# -ne 4 ]
then
	echo "Arguments missing"
	exit
fi

wd="$(pwd)"

#condition="$1"
#vcfs="$2"
#rep="$3"
#inv="$4"

## For testing
condition="Standard"
vcfs="Standard"
rep="Rep_6"
inv="HsInv0015"


pops="EUR AFR EAS"

mkdir -p $wd/Results_Minimac4/$condition/$rep/Inversions/
mkdir -p $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/

chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')

for pop in $pops
do

	## Create result file
	> $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/ImpGT_$pop

	## Sample Processing
	zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep "#CHROM" > $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/SamplesLine
	samples=$(cut -f10- $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/SamplesLine) 
	ncol=$(awk '{print NF}' $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/SamplesLine)
	impPos=$(zgrep "$inv" $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep -v "#" | cut -f2)

	for samp in $samples
	do
		## Prepare Reference panel
		bcftools view -s ^$samp --threads 10 -Ov -o $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Ref.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz

		cat $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Ref.vcf | \
		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Ref.vcf.gz
		tabix -f -p vcf $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Ref.vcf.gz
		
		minimac4 --compress-reference $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Ref.vcf.gz > $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Ref.msav
		

		## Prepare target
		bcftools view -s $samp --threads 10 -Ov -o $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Target.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz
		sed -i -E '/HsInv/ d' $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Target.vcf
		cat $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Target.vcf | \
    		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
     		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Target.vcf.gz
     		tabix -f -p vcf $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Target.vcf.gz
     		
     		
     		## Perform imputation
     		minimac4 \
     		--map $wd/../VCFs/30X/Maps/plink.chr${chr}.GRCh38.map \
     		--region "$chr:$(($impPos-1))-$(($impPos+1))" \
     		--threads 5 \
     		--output $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Out.vcf.gz \
     		--output-format vcf.gz \
     		--format GT,GP \
     		$wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Ref.msav \
     		$wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Target.vcf.gz


		

		echo "$samp" $(zcat $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Out.vcf.gz | grep "HsInv" | cut -f10 | sed 's/:/\t/g') | sed 's/ /\t/g' >> $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/ImpGT_$pop

#		rm $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Out.vcf.gz $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Out.log $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Ref* $wd/Results_Minimac4/$condition/$rep/Inversions/$inv/Target* 

	done
done


