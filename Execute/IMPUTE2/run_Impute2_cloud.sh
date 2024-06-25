#!/bin/bash

#SBATCH -p normal   # Partition to submit to
#SBATCH -J SplitImpute
#SBATCH -o /users/genomics/adriam/Ex_Def/Execute/Results_Impute2/logs/Split_%j_%A_%a.out 		# Sdt out file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID
#SBATCH -e /users/genomics/adriam/Ex_Def/Execute/Results_Impute2/logs/Split_%j_%A_%a.err 		# Sdt err file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID

wd="$(pwd)"


pops="EUR AFR EAS"
condition="$1"
vcfs="$2"
rep="$3"
invs=$(ls $wd/../VCFs/Inversion_Regions/$vcfs)
#pops="EUR"
#condition="Standard"
#vcfs="Standard"
#rep="Rep_1"
#invs="HsInv0124"

mkdir -p $wd/Results_Impute2
mkdir -p $wd/Results_Impute2/logs

> $wd/Results_Impute2/SampComp

for inv in $invs
do
	for pop in $pops
	do
		cat $wd/../VCFs/Inversion_Regions/$vcfs/$inv/${pop}_commSamps | tr " " "\n" | awk -v inv=$inv -v pop=$pop '{print inv"\t"pop"\t"$0}' >> $wd/Results_Impute2/SampComp
	done
done

nsamp=$(cat $wd/Results_Impute2/SampComp | wc -l)

# Número total de trabajos
total_jobs=$nsamp

rm -rf $wd/Results_Impute2/Chunks

mkdir -p $wd/Results_Impute2/Chunks

split -l 800 $wd/Results_Impute2/SampComp $wd/Results_Impute2/Chunks/chunk_
echo $total_jobs

for ch in $(ls $wd/Results_Impute2/Chunks)
do
  nch=$(cat $wd/Results_Impute2/Chunks/$ch | wc -l)
  sbatch --array=1-$nch $wd/IMPUTE2/run_Impute2_clouded.sh $wd/Results_Impute2/Chunks/$ch $condition $vcfs $rep

done
