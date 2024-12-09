#!/bin/bash

cd "${SLURM_SUBMIT_DIR}"

machine="${1}"
tool="${2}"
ecf_file="${3}"
output_path="${4}"
uniform_clusters_domains="${5}"
element_type="${6}"
domains_x="${7}"
domains_y="${8}"
domains_z="${9}"
elements_x="${10}"
elements_y="${11}"
elements_z="${12}"
dual_operator="${13}"

shift 3





mkdir -p "${output_path}"

infofile="${output_path}/info.txt"
echo -n > "${infofile}"
echo -n "datestr " >> "${infofile}"; date +%Y%m%d_%H%M%S >> "${infofile}"
echo "machine ${machine}" >> "${infofile}"
echo "tool ${tool}" >> "${infofile}"
echo -n "node " >> "${infofile}"; hostname >> "${infofile}"

if [[ "${dual_operator}" == *"GPU"* ]]; then
    if [ "${machine}" == "karolina" ] || [ "${machine}" == "e4red" ]; then
        echo -n "gpu " >> "${infofile}"; nvidia-smi -L >> "${infofile}"
    elif [ "${machine}" == "lumi" ]; then
        gpuinfo=$(rocm-smi --showproductname | grep GPU | head -n 1 | cut -d: -f3);
        echo -n "gpu " >> "${infofile}"; echo $gpuinfo >> "${infofile}"
    fi
else
    echo "cpu" >> "${infofile}"
fi

echo -n "cpus " >> "${infofile}"; numactl -s | grep physcpubind | cut -d' ' -f2- >> "${infofile}"
echo "ecf_file ${ecf_file}" >> "${infofile}"
echo "output_path ${output_path}" >> "${infofile}"
echo "uniform_clusters_domains ${uniform_clusters_domains}" >> "${infofile}"
echo "element_type ${element_type}" >> "${infofile}"
echo "domains_x ${domains_x}" >> "${infofile}"
echo "domains_y ${domains_y}" >> "${infofile}"
echo "domains_z ${domains_z}" >> "${infofile}"
echo "elements_x ${elements_x}" >> "${infofile}"
echo "elements_y ${elements_y}" >> "${infofile}"
echo "elements_z ${elements_z}" >> "${infofile}"
echo "dual_operator ${dual_operator}" >> "${infofile}"





command=""
if [ "${machine}" == "karolina" ]; then
    if [ "${tool}" == "cudalegacy" ]; then
        source env/it4i.karolina.cuda.gcc.32.sh legacy
        command="mpirun -n 1 --bind-to numa ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    elif [ "${tool}" == "cudamodern" ]; then
        source env/it4i.karolina.cuda.gcc.32.sh modern
        command="mpirun -n 1 --bind-to numa ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    elif [ "${tool}" == "suitesparse" ]; then
        source env/it4i.karolina.cuda.gcc.32.sh legacy
        command="mpirun -n 1 --bind-to numa ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    elif [ "${tool}" == "mklpardiso" ]; then
        source env/it4i.karolina.intel.32.sh
        export OMP_NUM_THREADS="$(nproc),1"
        command="mpirun -n 1 ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    else
        echo wrong tool
        exit 72
    fi
elif [ "${machine}" == "lumi" ]; then
    if [ "${tool}" == "rocm" ]; then
        source env/csc.lumi.rocm.mpich.sh
        command="srun -n 1 ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    elif [ "${tool}" == "suitesparse" ]; then
        source env/csc.lumi.rocm.mpich.sh
        command="srun -n 1 ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    elif [ "${tool}" == "mklpardiso" ]; then
        source env/csc.lumi.intel.sh
        export OMP_NUM_THREADS=7,1
        command="srun -n 1 ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    else
        echo wrong tool
        exit 72
    fi
elif [ "${machine}" == "e4red" ]; then
    if [ "${tool}" == "cudamodern" ]; then
        source env/e4.red.cuda.32.sh
        command="mpirun -n 1 --bind-to numa ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    elif [ "${tool}" == "suitesparse" ]; then
        source env/e4.red.cuda.32.sh
        export OMP_NUM_THREADS=72,1 # for cpu-only work, use all the threads
        command="mpirun -n 1 --bind-to numa ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    else
        echo wrong tool
        exit 72
    fi
elif [ "${machine}" == "tiber" ]; then
    if [ "${tool}" == "oneapi" ]; then
        source env/intel.tiber.oneapi.sh suitesparse
        command="mpirun -n 1 ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    elif [ "${tool}" == "suitesparse" ]; then
        source env/intel.tiber.oneapi.sh suitesparse
        command="mpirun -n 1 ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    elif [ "${tool}" == "mklpardiso" ]; then
        source env/intel.tiber.oneapi.sh mklpardiso
        command="mpirun -n 1 ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    else
        echo wrong tool
        exit 72
    fi
elif [[ "${machine}" == "spr"* ]]; then
    if [ "${tool}" == "suitesparse" ]; then
        source env/it4i.cs.spr.sh suitesparse
        command="./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    elif [ "${tool}" == "mklpardiso" ]; then
        source env/it4i.cs.spr.sh mklpardiso
        command="./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
    else
        echo wrong tool
        exit 72
    fi
else
    echo wrong machine
    exit 73
fi





timeout -v 300s bash -c "${command}" 2> "${output_path}/timeout.txt"
exitcode=$?





rm -f core*

exit ${exitcode}
