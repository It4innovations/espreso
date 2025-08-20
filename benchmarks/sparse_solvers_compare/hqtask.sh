#!/bin/bash

# assuming Karolina



taskdir="${HQ_ENTRY}"

if [ -f "${taskdir}/.finished" ]
then
    exit 0
fi

touch "${taskdir}/.started"

source "${taskdir}/config.sh"

ecf_file="${ESPRESO_BENCHMARK_ecf_file}"
ecf_file_name=$(basename "${ecf_file}" .ecf)

(
    echo -n "datestr: "; date +%Y%m%d_%H%M%S
    echo -n "hostname: "; hostname
    echo -n "numactl_cpus: "; numactl -s | grep physcpubind | cut -d' ' -f2-
    echo -n "nvidia_gpu: "; nvidia-smi -L
) > "${taskdir}/info.txt"

toexecute="mpirun -n 1 --bind-to numa ./build/espreso -c \"${ecf_file}\" >\"${taskdir}/stdout.txt\" 2>\"${taskdir}/stderr.txt\""

(
    cd "${taskdir}"
    ln -s "last/${ecf_file_name}.log" "log_out.txt"
    ln -s "last/${ecf_file_name}.err" "log_err.txt"
)



timeout -v 1800s bash -c "${toexecute}" 2> "${taskdir}/timeout.txt"
exitcode=$?



rm -f core*

echo "${exitcode}" > "${taskdir}/.finished"

exit ${exitcode}
