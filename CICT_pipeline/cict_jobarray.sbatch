#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=1800GB
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/as15096/eric/outputs/cict_par/cict_scaling/out_%A_%a.txt 
#SBATCH --error=/scratch/as15096/eric/outputs/cict_par/cict_scaling/err_%A_%a.txt

module purge
module load r/gcc/4.1.2

echo ${experiment_dir}
cd ${experiment_dir}


#Rscript /scratch/as15096/eric/Algorithms/CICT/runCICTEval2.R runCICT_par /scratch/as15096/eric/outputs/cict_par/cict_scaling/parConf_$SLURM_ARRAY_TASK_ID.yaml TRUE
command time -v -o /scratch/as15096/eric/outputs/cict_par/cict_scaling/time_$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID.txt Rscript /scratch/as15096/eric/Algorithms/CICT/runCICTEval2.R ${operation} ${experiment_dir}/parConf_$SLURM_ARRAY_TASK_ID.yaml TRUE


#command time -v -o time_$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID.txt Rscript /scratch/as15096/eric/Algorithms/CICT/runCICTEval2.R runCICT_par ${experiment_dir}/parConf_$SLURM_ARRAY_TASK_ID.yaml TRUE
