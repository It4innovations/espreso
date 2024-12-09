#!/bin/bash
#SBATCH --account OPEN-28-64
#SBATCH --partition p10-intel
#SBATCH --nodes 1
#SBATCH --time 1:00:00



workers_outdir_root="${SLURM_SUBMIT_DIR}/benchmarks/dualop_cpu_vs_gpu/workers_outerr"

hq_bin_dir="${SLURM_SUBMIT_DIR}/dependencies/HyperQueue"
HQBIN="${hq_bin_dir}/hq"

for i in $(seq -w 0 7)
do
    ii=$(( 10#${i} ))
    startcore=$(( ${ii}*12 ))
    endcore=$(( (${ii}+1)*12-1 ))
    outfile="${workers_outdir_root}/last/worker${i}_out.txt"
    errfile="${workers_outdir_root}/last/worker${i}_err.txt"
    numanode=$(( 8+${ii} ))
    numactl -C "${startcore}-${endcore}" -p "${numanode}" "${HQBIN}" worker start --idle-timeout=1m --heartbeat=5m > "${outfile}" 2> "${errfile}" &
    sleep 1
done

wait
