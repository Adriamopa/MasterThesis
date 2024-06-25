#!/bin/bash

if [ $# -ne 4 ]
then
	echo "Arguments missing"
	exit
fi

wd="$(pwd)"

condition="$1"
vcfs="$2"
rep="$3"
inv="$4"

## For testing
#condition="Standard"
#vcfs="Standard"
#inv="HsInv0015"


pops="EUR AFR EAS"

mkdir -p $wd/Results_Beagle/$condition/$rep/Inversions/
mkdir -p $wd/Results_Beagle/$condition/$rep/Inversions/$inv/

chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')


for pop in $pops
do

	## Create result file
	> $wd/Results_Beagle/$condition/$rep/Inversions/$inv/ImpGT_$pop

	## Sample Processing
	zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep "#CHROM" > $wd/Results_Beagle/$condition/$rep/Inversions/$inv/SamplesLine
	samples=$(cut -f10- $wd/Results_Beagle/$condition/$rep/Inversions/$inv/SamplesLine) 
	ncol=$(awk '{print NF}' $wd/Results_Beagle/$condition/$rep/Inversions/$inv/SamplesLine)

	for samp in $samples
	do

		## Prepare Reference panel
		bcftools view -s ^$samp --threads 10 -Ov -o $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Ref.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz

		cat $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Ref.vcf | \
		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Ref.vcf.gz
		tabix -f -p vcf $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Ref.vcf.gz
		

		## Prepare target
		bcftools view -s $samp --threads 10 -Ov -o $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Target.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz
		sed -i -E '/HsInv/ d' $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Target.vcf
		cat $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Target.vcf | \
    		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
     		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Target.vcf.gz
     		tabix -f -p vcf $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Target.vcf.gz

		if [ $chr == "X" ]
		then
			java -Xmx3g -jar $wd/../Programes/imputation_methods/beagle.22Jul22.46e.jar \
			ref=$wd/Results_Beagle/$condition/$rep/Inversions/$inv/Ref.vcf \
			gt=$wd/Results_Beagle/$condition/$rep/Inversions/$inv/Target.vcf \
			out=$wd/Results_Beagle/$condition/$rep/Inversions/$inv/Out \
			map=$wd/../VCFs/30X/Maps/plink.chr${chr}.GRCh38.map \
			ne=15000 \
			gp=true \
			em=false
		else
			java -Xmx3g -jar $wd/../Programes/imputation_methods/beagle.22Jul22.46e.jar \
			ref=$wd/Results_Beagle/$condition/$rep/Inversions/$inv/Ref.vcf \
			gt=$wd/Results_Beagle/$condition/$rep/Inversions/$inv/Target.vcf \
			out=$wd/Results_Beagle/$condition/$rep/Inversions/$inv/Out \
			map=$wd/../VCFs/30X/Maps/plink.chr${chr}.GRCh38.map \
			ne=20000 \
			gp=true \
			em=false
		fi

		# java -Xmx3g -jar $wd/../Software/beagle.05May22.33a.jar \
		# 	ref=$wd/Results/$inv/Ref.vcf \
		# 	gt=$wd/Results/$inv/Target.vcf \
		# 	out=$wd/Results/$inv/Out \
		# 	map=$wd/Maps/plink.chr${chr}.GRCh37.map

		echo "$samp" $(zcat $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Out.vcf.gz | grep "HsInv" | cut -f10 | sed 's/:/\t/g') | sed 's/ /\t/g' >> $wd/Results_Beagle/$condition/$rep/Inversions/$inv/ImpGT_$pop

		rm $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Out.vcf.gz $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Out.log $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Ref* $wd/Results_Beagle/$condition/$rep/Inversions/$inv/Target*

	done
done


