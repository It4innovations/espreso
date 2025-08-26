#!/bin/bash



taskid=0

function create_task {

    tasknum=$(printf "%06d" ${taskid})
    taskdir="${tasksdir}/${tasknum}"
    mkdir -p "${taskdir}"
    configfile="${taskdir}/config.sh"
    taskid=$((taskid+1))

    (
        echo "#!/bin/bash"
        echo
        echo "source env/it4i.karolina.gcc.cuda.mkl.ss.sh legacy"
        echo
        echo "export ESPRESO_BENCHMARK_ecf_file=\"${ecf_file}\""
        echo
        echo "export ECF_OUTPUT_PATH=${taskdir}"
        echo "export ECF_NUM_DOMAINS=${num_domains}"
        echo "export ECF_GPU_TRSM_TRIA_STRATEGY=${splitting_strategy}"
        echo "export ECF_PARTITION_STRATEGY=${partition_strategy}"
        echo "export ECF_CHUNK_COUNT=${chunk_count}"
    ) > "${configfile}"
    chmod +x "${configfile}"
}





basedir="benchmarks/dualop_cuda_spdp_real_problems"

if [ ! -d "${basedir}" ]
then
    echo "must be run from espreso root directory"
    exit 7
fi

datestr="$(date +%Y%m%d_%H%M%S)"

rundir="${basedir}/runs/${datestr}"
mkdir -p "${rundir}"

tasksdir="${rundir}/tasks"
mkdir -p "${tasksdir}"






array_ecf_files=("${basedir}/cooler.ecf" "${basedir}/drillbit.ecf" "${basedir}/wheel.ecf")




settings_tuples=(
    "AUTO      AUTO        0"
    "SPLIT_RHS CHUNK_COUNT 1"
)

num_domains_array=(16 32 64 128 256 512 1024 2048 4096 8192)



for ecf_file in "${array_ecf_files[@]}"
do
    for tuple in "${settings_tuples[@]}"
    do
        set -- $tuple
        splitting_strategy=$1
        partition_strategy=$2
        chunk_count=$3

        for num_domains in "${num_domains_array[@]}"
        do

            create_task

        done

    done
done



echo "created ${taskid} tasks"

echo "submit tasks to HQ using:"
echo "    ./benchmarks/dualop_cuda_spdp_real_problems/submit_tasks.sh ${rundir}"
