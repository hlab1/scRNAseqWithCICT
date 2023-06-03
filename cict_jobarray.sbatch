#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50GB
#SBATCH --time=10:00
#SBATCH --output=/scratch/as15096/eric/outputs_v2/cict_par/sens_modeling_choices/out_%A_%a.txt 
#SBATCH --error=/scratch/as15096/eric/outputs_v2/cict_par/sens_modeling_choices/err_%A_%a.txt

module purge
module load r/gcc/4.1.2

echo ${experiment_dir}
cd ${experiment_dir}


#Rscript /scratch/as15096/eric/Algorithms/CICT/runCICTEval2.R runCICT_par /scratch/as15096/eric/outputs_v2/cict_par/sens_modeling_choices/parConf_$SLURM_ARRAY_TASK_ID.yaml TRUE

command time -v -o /scratch/as15096/eric/outputs_v2/cict_par/sens_modeling_choices/time_$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID.txt Rscript /scratch/as15096/eric/Algorithms/CICT/runCICTEval2.R ${operation} ${experiment_dir}/parConf_$SLURM_ARRAY_TASK_ID.yaml TRUE


#command time -v -o time_$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID.txt Rscript /scratch/as15096/eric/Algorithms/CICT/runCICTEval2.R runCICT_par ${experiment_dir}/parConf_$SLURM_ARRAY_TASK_ID.yaml TRUE
