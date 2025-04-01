#!/bin/bash
#SBATCH --account=ehpc202
#SBATCH --qos=acc_ehpc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00

hqbin="dependencies/HyperQueue/hq"

export ROCR_VISIBLE_DEVICES=""
numactl -C +0-19 ${hqbin} worker start --idle-timeout=1m
