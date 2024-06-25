#!/bin/bash

#SBATCH -J Full_Impute2
#SBATCH -o /users/genomics/adriam/Ex_Def/Execute/Results_Impute2/logs/Slurm_%j_%A_%a.out 		# Sdt out file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID
#SBATCH -e /users/genomics/adriam/Ex_Def/Execute/Results_Impute2/logs/Slurm_%j_%A_%a.err 		# Sdt err file name: %j SLURM_JOB_ID; %A SLURM_ARRAY_JOB_ID, %a SLURM_ARRAY_TASK_ID

wd=$(pwd)

narr=$(ls $wd/../VCFs/Inversion_Regions/Standard/ | wc -l)

run_id=$(sbatch --array=1-${narr} --parsable $wd/IMPUTE2/primer.sh)

sbatch --dependency=afterok:${run_id} $wd/IMPUTE2/segon.sh
