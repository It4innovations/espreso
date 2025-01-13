#!/bin/bash
#SBATCH --account eupuser
#SBATCH --partition t-gpu-gh200
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=72
#SBATCH --exclusive
#SBATCH --time 12:00:00



if nvidia-smi | grep "ERR!" >/dev/null 2>/dev/null
then
    echo "weird gpu, ending this job"
    nvidia-smi
    exit 1
fi



export ROCR_VISIBLE_DEVICES=""

hq_bin_dir="${SLURM_SUBMIT_DIR}/dependencies/HyperQueue"
HQBIN="${hq_bin_dir}/hq"

${HQBIN} worker start --idle-timeout=1m
