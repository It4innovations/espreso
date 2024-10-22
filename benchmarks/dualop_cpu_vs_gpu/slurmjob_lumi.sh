#!/bin/bash
#SBATCH --account project_465001013
#SBATCH --partition small-g
#SBATCH --ntasks 1
#SBATCH --gpus-per-task 1
#SBATCH --cpus-per-task=7
#SBATCH --time 4:00:00



hq_bin_dir="${SLURM_SUBMIT_DIR}/dependencies/HyperQueue"
HQBIN="${hq_bin_dir}/hq"

numactl -C +0-6 "${HQBIN}" worker start --idle-timeout=1m
