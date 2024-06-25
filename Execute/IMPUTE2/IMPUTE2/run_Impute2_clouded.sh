#!/bin/bash

#SBATCH -p normal   # Partition to submit to
#SBATCH --cpus-per-task=1	#number-cpus per task
#SBATCH --mem-per-cpu 8Gb      # Memory in MB
#SBATCH -J Impute             # job name
#SBATCH -o /users/genomics/adriam/Ex_Def/Execute/Results_Impute2/logs/Imp_%j_%A_%a.out 		# Sdt out file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID
#SBATCH -e /users/genomics/adriam/Ex_Def/Execute/Results_Impute2/logs/Imp_%j_%A_%a.err 		# Sdt err file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID

wd=$(pwd)

module load plink/1.9
module load R/4.3.3-foss-2020b
module load BCFtools/1.12-GCC-10.2.0
module load tabix/0.2.6-GCCcore-10.2.0
module load vcftools/0.1.16
module load Perl/5.32.0-GCCcore-10.2.0

chunk_f=$1
condition=$2
vcfs=$3
rep=$4

declare -a invs
declare -a pops
declare -a samps

while read col1 col2 col3
do
	invs=(${invs[@]} $col1)
	pops=(${pops[@]} $col2)
	samps=(${samps[@]} $col3)

done < $chunk_f

i=$(($SLURM_ARRAY_TASK_ID -1))

inv=${invs[i]}
pop=${pops[i]}
samp=${samps[i]}

rdir="$wd/Results_Impute2/$condition/$rep/Inversions/$inv/$pop/$samp"
mkdir -p $rdir

vcf="$wd/../VCFs/Inversion_Regions/$vcfs/$inv/$pop.vcf.gz"

zcat $vcf | grep "#CHROM" > ${rdir}/SamplesLine

chr=$(grep "$inv," $wd/../VCFs/30X/InvCoordenates_hg38_v1.csv | cut -d',' -f4 | sed 's/chr//g')
impPos=$(zgrep "$inv" $vcf | grep -v "#" | cut -f2)

## Prepare Reference
bcftools view -s ^$samp --threads 10 -Ov -o ${rdir}/Ref.vcf $vcf

cat ${rdir}/Ref.vcf | sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g'  | sed -e 's/\t0$/\t0|0/g' | bgzip > ${rdir}/Ref.vcf.gz
tabix -f -p vcf ${rdir}/Ref.vcf.gz

$wd/../Programes/vcf2impute_gen.pl -vcf ${rdir}/Ref.vcf -gen ${rdir}/Ref.gen


## Prepare Target
bcftools view -s $samp --threads 10 $vcf | grep -v $inv > ${rdir}/Target.vcf

cat ${rdir}/Target.vcf | sed 's/\t0\t/\t0|0\t/g' | sed -e 's/\t0\t/\t0|0\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1\t/\t1|1\t/g' | sed -e 's/\t1$/\t1|1/g' | sed -e 's/\t0$/\t0|0/g' | bgzip > ${rdir}/Target.vcf.gz
tabix -f -p vcf ${rdir}/Target.vcf.gz

$wd/../Programes/vcf2impute_gen.pl -vcf ${rdir}/Target.vcf -gen ${rdir}/Target.gen


## Run IMPUTE2
if [ $chr == "X" ]
then
	Rscript $wd/IMPUTE2/MakeForX_Impute2.R $inv $pop $samp $condition $rep

	"$wd/../Programes/imputation_methods/impute_v2.3.2_x86_64_dynamic/impute2" \
	-g_ref ${rdir}/Ref.gen.gz \
	-g ${rdir}/Target.gen.gz \
	-int $impPos $impPos \
	-buffer 200 kb \
	-o ${rdir}/Out \
	-m $wd/../VCFs/30X/Maps/chrX_hg38_genmap.txt \
	-sample_g_ref ${rdir}/Ref.sample \
	-sample_g ${rdir}/Target.sample \
	-Ne 20000 \
	-chrX
rm ${rdir}/Ref.sample ${rdir}/Target.sample

else
	"$wd/../Programes/imputation_methods/impute_v2.3.2_x86_64_dynamic/impute2" \
	-g_ref ${rdir}/Ref.gen.gz \
	-g ${rdir}/Target.gen.gz \
	-int $impPos $impPos \
	-buffer 200 kb \
	-o ${rdir}/Out \
	-m $wd/../VCFs/30X/Maps/chr${chr}_hg38_genmap.txt  \
	-Ne 20000
fi

## Save the results
echo "$samp" $(grep "HsInv" ${rdir}/Out | cut -d" " -f6- | sed 's/ /\t/g') | sed 's/ /\t/g' >> $wd/Results_Impute2/$condition/$rep/Inversions/$inv/ImpGT_$pop

rm ${rdir}/Ref* ${rdir}/Target*
rm $wd/Results_Impute2/logs/*
