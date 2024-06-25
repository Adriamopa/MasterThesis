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
#inv="HsInv0015"
#condition="Standard"
#vcfs="Standard"
#rep="rep_1"

chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')

BP1=$(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f5)-200000))
BP2=$(($(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f8)+200000))

pops="EUR AFR EAS"

mkdir -p $wd/Results_Impute5/$condition/$rep/Inversions
mkdir -p $wd/Results_Impute5/$condition/$rep/Inversions/$inv

if [ $chr == "X" ]
then
	for pop in $pops
	do

		## Create result file
		> $wd/Results_Impute5/$condition/$rep/Inversions/$inv/ImpGT_$pop

		## Sample Processing
		zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep "#CHROM" > $wd/Results_Impute5/$condition/$rep/Inversions/$inv/SamplesLine
		samples=$(cut -f10- $wd/Results_Impute5/$condition/$rep/Inversions/$inv/SamplesLine) 
		echo "$samples" | sed 's/\t/\n/g' > $wd/Results_Impute5/$condition/$rep/Inversions/$inv/indivs.txt
		grep -F -f $wd/Results_Impute5/$condition/$rep/Inversions/$inv/indivs.txt $wd/../VCFs/30X/Panel30x | cut -f1,4 > $wd/Results_Impute5/$condition/$rep/Inversions/$inv/SexedSamples
		ncol=$(awk '{print NF}' $wd/Results_Impute5/$condition/$rep/Inversions/$inv/SamplesLine)
		impPos=$(zgrep "$inv" $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep -v "#" | cut -f2)
		
		#samples="NA12878"

		for samp in $samples
		do
			sex="$(grep "$samp" $wd/Results_Impute5/$condition/$rep/Inversions/$inv/SexedSamples | cut -f2)"


			## Prepare Reference panel
			bcftools view -s ^$samp --threads 10 -Ov -o $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz

			cat $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf | \
    		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
     		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz
			tabix -f -p vcf $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz
			

			## Prepare target
			bcftools view -s $samp --threads 10 -Ov -o $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz
			sed -i -E '/HsInv/ d' $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf
			cat $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf | \
    		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
     		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz
	     		tabix -f -p vcf $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz

			if [ $sex = "male" ]
			then
			impute5 \
				--h $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz \
				--m $wd/../VCFs/30X/Maps/chrX.b38.gmap.gz \
				--g $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz \
				--r X:$(($impPos-1))-$(($impPos+1)) \
				--o $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Out.vcf.gz \
				--l $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Log \
				--pbwt-cm 0.005 \
				--threads 10 \
				--ne 20000 \
				--buffer-region X:$BP1-$BP2 \
				--no-out-index \
				--haploid
				
			else
			impute5 \
				--h $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz \
				--m $wd/../VCFs/30X/Maps/chrX.b38.gmap.gz \
				--g $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz \
				--r X:$(($impPos-1))-$(($impPos+1)) \
				--o $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Out.vcf.gz \
				--l $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Log \
				--pbwt-cm 0.005 \
				--threads 10 \
				--ne 20000 \
				--buffer-region X:$BP1-$BP2 \
				--no-out-index
			fi

			echo "$samp" $(zcat $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Out.vcf.gz | grep "HsInv" | cut -f10 | sed 's/:/\t/g') | sed 's/ /\t/g' >> $wd/Results_Impute5/$condition/$rep/Inversions/$inv/ImpGT_$pop

			rm $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Out.vcf.gz $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz.tbi \
			  $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz.tbi $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Header

		done
	done

else

	for pop in $pops
	do

		## Create result file
		> $wd/Results_Impute5/$condition/$rep/Inversions/$inv/ImpGT_$pop

		## Sample Processing
		zcat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep "#CHROM" > $wd/Results_Impute5/$condition/$rep/Inversions/$inv/SamplesLine
		samples=$(cut -f10- $wd/Results_Impute5/$condition/$rep/Inversions/$inv/SamplesLine) 
		ncol=$(awk '{print NF}' $wd/Results_Impute5/$condition/$rep/Inversions/$inv/SamplesLine)
		impPos=$(zgrep "$inv" $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz | grep -v "#" | cut -f2)
		
		#samples="NA12878"

		for samp in $samples
		do
		
			## Prepare Reference panel
			bcftools view -s ^$samp --threads 10 -Ov -o $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz

	    		cat $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf | \
    		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
     		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz
	    		tabix -f -p vcf $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz
			

			## Prepare target
			bcftools view -s $samp --threads 10 -Ov -o $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf $wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz
			sed -i -E '/HsInv/ d' $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf
			cat $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf | \
    		sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | \
     		sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz
			tabix -f -p vcf $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz

			impute5 \
				--h $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz \
				--m $wd/../VCFs/30X/Maps/chr${chr}.b38.gmap.gz \
				--g $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz \
				--r $chr:$(($impPos-1))-$(($impPos+1)) \
				--o $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Out.vcf.gz \
				--l $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Log \
				--pbwt-cm 0.005 \
				--threads 10 \
				--ne 20000 \
				--buffer-region $chr:$BP1-$BP2 \
				--no-out-index

			echo "$samp" $(zcat $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Out.vcf.gz | grep "HsInv" | cut -f10 | sed 's/:/\t/g') | sed 's/ /\t/g' >> $wd/Results_Impute5/$condition/$rep/Inversions/$inv/ImpGT_$pop

			rm $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Out.vcf.gz $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Ref.vcf.gz.tbi \
			  $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz $wd/Results_Impute5/$condition/$rep/Inversions/$inv/Target.vcf.gz.tbi
		done
	done


fi

