#!/bin/bash
#SBATCH --account OPEN-28-64
#SBATCH --partition qgpu
#SBATCH --nodes 1
#SBATCH --gpus-per-node=1
#SBATCH --time 4:00:00



if nvidia-smi | grep "ERR!" >/dev/null 2>/dev/null
then
    echo "weird gpu, ending this job"
    nvidia-smi
    exit 1
fi



module load HyperQueue/0.19.0

export ROCR_VISIBLE_DEVICES=""

hq worker start --idle-timeout=1m
