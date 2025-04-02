#!/bin/bash

taskdir="${HQ_ENTRY}"

if [ -f "${taskdir}/.finished" ]
then
    exit 0
fi

touch "${taskdir}/.started"

source "${taskdir}/config.sh"

machine="${ESPRESO_BENCHMARK_machine}"
ecf_file="${ESPRESO_BENCHMARK_ecf_file}"
ecf_file_name=$(basename "${ecf_file}" .ecf)

(
    echo -n "datestr: "; date +%Y%m%d_%H%M%S
    echo -n "hostname: "; hostname
    echo -n "numactl_cpus: "; numactl -s | grep physcpubind | cut -d' ' -f2-
    if [ "${machine}" == "karolina" ]; then
        echo -n "nvidia_gpu: "; nvidia-smi -L
    fi
    if [ "${machine}" == "mn5" ]; then
        echo -n "nvidia_gpu: "; nvidia-smi -L
    fi
) > "${taskdir}/info.txt"

ecf_args="${taskdir} TRUE ${ESPRESO_BENCHMARK_element_type} ${ESPRESO_BENCHMARK_domains_x} ${ESPRESO_BENCHMARK_domains_y} ${ESPRESO_BENCHMARK_domains_z} ${ESPRESO_BENCHMARK_elements_x} ${ESPRESO_BENCHMARK_elements_y} ${ESPRESO_BENCHMARK_elements_z} ${ESPRESO_BENCHMARK_dual_operator}"

toexecute=""
if [ "${machine}" == "karolina" ]; then
    toexecute="mpirun -n 1 --bind-to numa ./build/espreso -c \"${ecf_file}\" ${ecf_args} >\"${taskdir}/stdout.txt\" 2>\"${taskdir}/stderr.txt\""
fi
if [ "${machine}" == "mn5" ]; then
    toexecute="mpirun -n 1 --bind-to numa ./build/espreso -c \"${ecf_file}\" ${ecf_args} >\"${taskdir}/stdout.txt\" 2>\"${taskdir}/stderr.txt\""
fi

(
    cd "${taskdir}"
    ln -s "last/${ecf_file_name}.log" "out.txt"
    ln -s "last/${ecf_file_name}.err" "err.txt"
)



timeout -v 300s bash -c "${toexecute}" 2> "${taskdir}/timeout.txt"
exitcode=$?



rm -f core*

echo "${exitcode}" > "${taskdir}/.finished"

exit ${exitcode}
