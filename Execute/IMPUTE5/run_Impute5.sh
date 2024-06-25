#!/bin/bash

## DON'T RUN, run  run_combined_Impute5.sh instead

if [ $# -ne 3 ]
then
	echo "Arguments missing"
	exit
fi

wd="$(pwd)"

condition="$1"
vcfs="$2"
inv="$3"

## For testing
#inv="HsInv0015"
#condition="Standard"
#vcfs="Standard"

chr=$(grep "$inv[[:space:]]" $wd/../VCFs/30X/Common_hg38 | cut -f2 | sed 's/chr//g')

if [[ $chr =~ [Xx] ]]
then
   bash run_chrX_Impute5.sh $condition $vcfs $inv
   exit
fi

BP1=$(($(grep "$inv[[:space:]]" $wd/../VCFs/30X/Common_hg38 | cut -f3)-200000))
BP2=$(($(grep "$inv[[:space:]]" $wd/../VCFs/30X/Common_hg38 | cut -f6)+200000))

pops="EUR AFR EAS"

mkdir -p $wd/Results_Impute5/Inversions/$condition
mkdir -p $wd/Results_Impute5/Inversions/$condition/$inv

for pop in $pops
do

	## Create result file
	> $wd/Results_Impute5/Inversions/$condition/$inv/ImpGT_$pop

	## Sample Processing
	zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep "#CHROM" > $wd/Results_Impute5/Inversions/$condition/$inv/SamplesLine
	samples=$(cut -f10- $wd/Results_Impute5/Inversions/$condition/$inv/SamplesLine) 
	ncol=$(awk '{print NF}' $wd/Results_Impute5/Inversions/$condition/$inv/SamplesLine)
	impPos=$(zgrep "$inv" $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep -v "#" | cut -f2)
	
	#samples="NA12878"

	for samp in $samples
	do
		selcol="$(sed 's/\t/\n/g' $wd/Results_Impute5/Inversions/$condition/$inv/SamplesLine | cat -n | grep $samp | awk '{print $1}')"

    		zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep "##" > $wd/Results_Impute5/Inversions/$condition/$inv/Header

		## Prepare ref
		if [ $selcol == $ncol ]
		then
			zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep -v "##" | cut -f1-$(($selcol-1)) > $wd/Results_Impute5/Inversions/$condition/$inv/Ref.vcf
		else
			zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep -v "##" | cut -f1-$(($selcol-1)),$(($selcol+1))- > $wd/Results_Impute5/Inversions/$condition/$inv/Ref.vcf
		fi

    		cat $wd/Results_Impute5/Inversions/$condition/$inv/Header $wd/Results_Impute5/Inversions/$condition/$inv/Ref.vcf | \
      		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
      		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Impute5/Inversions/$condition/$inv/Ref.vcf.gz
    		tabix -f -p vcf $wd/Results_Impute5/Inversions/$condition/$inv/Ref.vcf.gz
		

		## Prepare target
		zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep -v "##" | cut -f1-9,$selcol | grep -v "HsInv" > $wd/Results_Impute5/Inversions/$condition/$inv/Target.vcf
		
		cat $wd/Results_Impute5/Inversions/$condition/$inv/Header $wd/Results_Impute5/Inversions/$condition/$inv/Target.vcf | \
      		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
      		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Impute5/Inversions/$condition/$inv/Target.vcf.gz
		tabix -f -p vcf $wd/Results_Impute5/Inversions/$condition/$inv/Target.vcf.gz

		impute5 \
			--h $wd/Results_Impute5/Inversions/$condition/$inv/Ref.vcf.gz \
			--m $wd/../VCFs/30X/Maps/chr${chr}.b38.gmap.gz \
			--g $wd/Results_Impute5/Inversions/$condition/$inv/Target.vcf.gz \
			--r $chr:$(($impPos-1))-$(($impPos+1)) \
			--o $wd/Results_Impute5/Inversions/$condition/$inv/Out.vcf.gz \
			--l $wd/Results_Impute5/Inversions/$condition/$inv/Log \
			--pbwt-cm 0.005 \
			--threads 10 \
			--ne 20000 \
			--buffer-region $chr:$BP1-$BP2 \
			--no-out-index

		echo "$samp" $(zcat $wd/Results_Impute5/Inversions/$condition/$inv/Out.vcf.gz | grep "HsInv" | cut -f10 | sed 's/:/\t/g') | sed 's/ /\t/g' >> $wd/Results_Impute5/Inversions/$condition/$inv/ImpGT_$pop

		rm $wd/Results_Impute5/Inversions/$condition/$inv/Out.vcf.gz $wd/Results_Impute5/Inversions/$condition/$inv/Ref.vcf $wd/Results_Impute5/Inversions/$condition/$inv/Ref.vcf.gz $wd/Results_Impute5/Inversions/$condition/$inv/Ref.vcf.gz.tbi \
		  $wd/Results_Impute5/Inversions/$condition/$inv/Target.vcf $wd/Results_Impute5/Inversions/$condition/$inv/Target.vcf.gz $wd/Results_Impute5/Inversions/$condition/$inv/Target.vcf.gz.tbi $wd/Results_Impute5/Inversions/$condition/$inv/Header

	done
done


