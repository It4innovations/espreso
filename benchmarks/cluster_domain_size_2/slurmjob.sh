#!/bin/bash
#SBATCH --account=ehpc202
#SBATCH --qos=acc_ehpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:4
#SBATCH --threads-per-core=1
#SBATCH --time=4:00:00

# assuming MN5



hqbin="dependencies/HyperQueue/hq"

export ROCR_VISIBLE_DEVICES=""

${hqbin} worker start --cpus=20 --idle-timeout=1m &
sleep 1
${hqbin} worker start --cpus=20 --idle-timeout=1m &
sleep 1
${hqbin} worker start --cpus=20 --idle-timeout=1m &
sleep 1
${hqbin} worker start --cpus=20 --idle-timeout=1m &
sleep 1

wait
