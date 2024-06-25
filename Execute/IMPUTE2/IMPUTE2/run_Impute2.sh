#!/bin/bash

if [ $# -ne 3 ]
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
#inv="HsInv0015"
#condition="Standard" or "Unphased"
#vcfs="Standard"
chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')

BP1=$(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f5)-200000))
BP2=$(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f8)+200000))

pops="EUR AFR EAS"

mkdir -p $wd/Results_Impute2
mkdir -p $wd/Results_Impute2/Inversions/$condition/$rep
mkdir -p $wd/Results_Impute2/Inversions/$condition/$rep/$inv

for pop in $pops
do

	## Create result file
	> $wd/Results_Impute2/Inversions/$condition/$rep/$inv/ImpGT_$pop

	## Sample Processing
	zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep "#CHROM" > $wd/Results_Impute2/Inversions/$condition/$rep/$inv/SamplesLine
	samples=$(cut -f10- $wd/Results_Impute2/Inversions/$condition/$rep/$inv/SamplesLine) 
	ncol=$(awk '{print NF}' $wd/Results_Impute2/Inversions/$condition/$rep/$inv/SamplesLine)

	impPos=$(zgrep "$inv" $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep -v "#" | cut -f2)

	echo $impPos

	for samp in $samples
	do

		## Prepare Reference panel
		bcftools view -s ^$samp --threads 10 -Ov -o $wd/Results_Impute2/$condition/$rep/$rep/Inversions/$inv/Ref.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz

		cat $wd/Results_Impute2/$condition/$rep/$rep/Inversions/$inv/Ref.vcf | \
		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Impute2/$condition/$rep/$rep/Inversions/$inv/Ref.vcf.gz
		tabix -f -p vcf $wd/Results_Impute2/$condition/$rep/$rep/Inversions/$inv/Ref.vcf.gz
		
		$wd/../Programes/vcf2impute_gen.pl -vcf $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.vcf -gen $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.gen


		## Prepare target
		bcftools view -s $samp --threads 10 -Ov -o $wd/Results_Impute2/$condition/$rep/$rep/Inversions/$inv/Target.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz
		sed -i -E '/HsInv/ d' $wd/Results_Impute2/$condition/$rep/$rep/Inversions/$inv/Target.vcf
		cat $wd/Results_Impute2/$condition/$rep/$rep/Inversions/$inv/Target.vcf | \
    		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
     		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Impute2/$condition/$rep/$rep/Inversions/$inv/Target.vcf.gz
     		tabix -f -p vcf $wd/Results_Impute2/$condition/$rep/$rep/Inversions/$inv/Target.vcf.gz

		$wd/../Programes/vcf2impute_gen.pl -vcf $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.vcf.gz -gen $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.gen


		if [ $chr == "X" ]
		then
			Rscript $wd/IMPUTE2/MakeForX_Impute2.R $inv $pop $samp
			
			impute2 \
				-g_ref $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.gen.gz \
				-g $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.gen.gz \
				-int $impPos $impPos \
				-buffer 200 kb \
				-o $wd/Results_Impute2/$condition/$rep/Inversions/$inv/Out \
				-m $wd/../VCFs/30X/Maps/chrX_hg38_genmap.txt \
				-sample_g_ref $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.sample \
				-sample_g $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.sample \
				-Ne 20000 \
				-chrX
			rm $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.sample $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.sample

		else
			impute2 \
				-g_ref $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.gen.gz \
				-g $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.gen.gz \
				-int $impPos $impPos \
				-buffer 200 kb \
				-o $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Out \
				-m $wd/../VCFs/30X/Maps/chr${chr}_hg38_genmap.txt  \
				-Ne 20000
		fi

		echo "$samp" $(grep "HsInv" $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Out | cut -d" " -f6- | sed 's/ /\t/g') | sed 's/ /\t/g' >> $wd/Results_Impute2/Inversions/$condition/$rep/$inv/ImpGT_$pop

		rm $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Out* $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.vcf $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.vcf.gz $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.vcf.gz.tbi \
		  $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.vcf $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.vcf.gz $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.vcf.gz.tbi \
		  $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Target.gen* $wd/Results_Impute2/Inversions/$condition/$rep/$inv/Ref.gen*

	done
done


