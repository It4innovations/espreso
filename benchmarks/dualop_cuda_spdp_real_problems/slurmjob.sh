#!/bin/bash
#SBATCH --account=OPEN-28-64
#SBATCH --partition=qgpu
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --time=4:00:00

# assuming Karolina



hqbin="dependencies/HyperQueue/hq"

export ROCR_VISIBLE_DEVICES=""
${hqbin} worker start --idle-timeout=1m
