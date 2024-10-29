#!/bin/bash
  
#SBATCH --account=phys033185
#SBATCH --job-name=submissionTest
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=0:0:10
#SBATCH --mem-per-cpu=100M 
#SBATCH --array=100-109

# Change to working directory, where the job was submitted from.
cd "${SLURM_SUBMIT_DIR}"

# Record some potentially useful details about the job:
#echo "Running on host $(hostname)"
#echo "Started on $(date)"
#echo "Directory is $(pwd)"
#echo "Slurm job ID is ${SLURM_JOBID}"
#echo "This jobs runs on the following machines:"
#echo "${SLURM_JOB_NODELIST}"
#printf "\n\n"

# Submit
#python ./programs/numpyLebwohlLasher.py 100 90 1 0
#python ./programs/cythonRunner.py 100 25 1 0
mpiexec -n 2 python ./programs/mpiLebwohlLasher.py 500 25 1 0

# Output the end time
#printf "\n\n"
#echo "Ended on: $(date)"

