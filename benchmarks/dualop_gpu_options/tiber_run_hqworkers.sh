#!/bin/bash

machine="tiber"

basedir="benchmarks/dualop_gpu_options"
if [ ! -f "${basedir}/${machine}_run_hqworkers.sh" ]
then
    echo "The script must be run from the espreso root folder"
    exit 1
fi

workers_outdir_root="${basedir}/workers_outerr"
if [[ ! -d "${workers_outdir_root}" ]]
then
    echo "Workers output directory '${workers_outdir_root}' does not exist"
    echo "Please create it"
    exit 2
fi

espreso_outdir_root="${basedir}/espreso_results"
if [[ ! -d "${espreso_outdir_root}" ]]
then
    echo "Espreso output directory '${espreso_outdir_root}' does not exist"
    echo "Please create it"
    exit 3
fi

hq_outdir_root="${basedir}/hq_outerr"
if [[ ! -d "${hq_outdir_root}" ]]
then
    echo "HyperQueue output directory '${hq_outdir_root}' does not exist"
    echo "Please create it"
    exit 4
fi

hq_bin_dir="dependencies/HyperQueue"
HQBIN="${hq_bin_dir}/hq"
if ! "${HQBIN}" server info > /dev/null 2> /dev/null
then
    echo "HyperQueue server is not running"
    echo "Start the server using 'numactl -C 0-5 ${HQBIN} server start' inside tmux"
    exit 5
fi







echo "starting workers"
echo

for i in $(seq -w 1 15) # 0 reserved for hq server and other stuff
do
    ii=$(( 10#${i} ))
    startcore=$(( ${ii}*6 ))
    endcore=$(( (${ii}+1)*6-1 ))
    export ONEAPI_DEVICE_SELECTOR="level_zero:${ii}"
    outfile="${workers_outdir_root}/last/worker${i}_out.txt"
    errfile="${workers_outdir_root}/last/worker${i}_err.txt"
    numactl -C "${startcore}-${endcore}" "${HQBIN}" worker start --time-limit=168h --idle-timeout=1m > "${outfile}" 2> "${errfile}" &
    sleep 1
done

echo
echo "workers started"
echo "waiting to finish"

wait

echo "workers finished"
