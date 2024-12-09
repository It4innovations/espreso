#!/bin/bash

machine="karolina"
# machine="lumi"
# machine="e4red"
# machine="tiber"
# machine="sprddr"
# machine="sprhbm"

# have to make sure the same version is compiled
# tool="cudalegacy"
# tool="cudamodern"
# tool="suitesparse"
tool="mklpardiso"
# tool="rocm"
# tool="oneapi"

dualoperators=""
# dualoperators="${dualoperators} EXPLICIT_GPU"
# dualoperators="${dualoperators} IMPLICIT_GPU"
dualoperators="${dualoperators} EXPLICIT_SC"
dualoperators="${dualoperators} IMPLICIT"



basedir="benchmarks/dualop_cpu_vs_gpu"
if [ ! -f "${basedir}/submit.sh" ]
then
    echo "The script must be run from the espreso root folder"
    exit 1
fi

if [ "${machine}" == "tiber" ] || [[ "${machine}" == "spr"* ]]; then
    workers_outdir_root="${basedir}/workers_outerr"
    if [[ ! -d "${workers_outdir_root}" ]]
    then
        echo "Workers output directory '${workers_outdir_root}' does not exist"
        echo "Please create it"
        exit 2
    fi
fi

if ! [ "${machine}" == "tiber" ]; then
    slurm_outdir_root="${basedir}/slurm_outerr"
    if [[ ! -d "${slurm_outdir_root}"   &&   ! ( -L "${slurm_outdir_root}" && -d "$(readlink ${slurm_outdir_root})" ) ]]
    then
        echo "Slurm output directory '${slurm_outdir_root}' does not exist"
        echo "Please create it"
        exit 2
    fi
fi

espreso_outdir_root="${basedir}/espreso_results"
if [[ ! -d "${espreso_outdir_root}"   &&   ! ( -L "${espreso_outdir_root}" && -d "$(readlink ${espreso_outdir_root})" ) ]]
then
    echo "Espreso output directory '${espreso_outdir_root}' does not exist"
    echo "Please create it"
    exit 3
fi

hq_outdir_root="${basedir}/hq_outerr"
if [[ ! -d "${hq_outdir_root}"   &&   ! ( -L "${hq_outdir_root}" && -d "$(readlink ${hq_outdir_root})" ) ]]
then
    echo "HyperQueue output directory '${hq_outdir_root}' does not exist"
    echo "Please create it"
    exit 4
fi

HQBIN=""
if [ "${machine}" == "karolina" ]; then
    module load HyperQueue/0.19.0
    HQBIN="hq"
else
    hq_bin_dir="dependencies/HyperQueue"
    HQBIN="${hq_bin_dir}/hq"

    # can download HQ using e.g.:
    # mkdir -p "${hq_bin_dir}"
    # cd "${hq_bin_dir}"
    # wget https://github.com/It4innovations/hyperqueue/releases/download/v0.19.0/hq-v0.19.0-linux-arm64-linux.tar.gz
    # tar -xf hq-v0.19.0-linux-arm64-linux.tar.gz
fi

if ! ${HQBIN} server info > /dev/null 2> /dev/null
then
    echo "HyperQueue server is not running"
    if [ "${machine}" == "tiber" ]; then
        echo "Start the server using 'numactl -C 0-5 ${HQBIN} server start' inside tmux"
    else
        echo "Start the server on this login node using 'hq server start' inside tmux"
    fi
    exit 5
fi



datestr="$(date +%Y%m%d_%H%M%S)"

hq_outdir="${hq_outdir_root}/${machine}_${tool}_${datestr}"
mkdir -p "${hq_outdir}"

espreso_outdir="${espreso_outdir_root}/${machine}_${tool}_${datestr}"
mkdir -p "${espreso_outdir}"

if [ "${machine}" == "tiber" ] || [[ "${machine}" == "spr"* ]]; then
    workers_outerr="${workers_outdir_root}/${machine}_${tool}_${datestr}"
    mkdir -p "${workers_outerr}"
    ln -f -s -T "${machine}_${tool}_${datestr}" "${workers_outdir_root}/last"
fi

if ! [ "${machine}" == "tiber" ]; then
    slurm_outdir="${slurm_outdir_root}/${machine}_${tool}_${datestr}"
    mkdir -p "${slurm_outdir}"
fi




num_cores_for_job=""
if [ "${machine}" == "karolina" ]; then
    num_cores_for_job="16"
elif [ "${machine}" == "lumi" ]; then
    num_cores_for_job="7"
elif [ "${machine}" == "e4red" ]; then
    num_cores_for_job="72"
elif [ "${machine}" == "tiber" ]; then
    num_cores_for_job="6"
elif [[ "${machine}" == "spr"* ]]; then
    num_cores_for_job="12"
fi


id_num=0
function submit_hq_job
{
    id_num=$((${id_num}+1))
    id_str=$(printf '%06d' ${id_num})
    output_path="${espreso_outdir}/${id_str}"
    echo "Submitting task ${id_str}"
    ${HQBIN} submit \
        --time-request=6m \
        --cpus="${num_cores_for_job} compact" \
        --stdout="${hq_outdir}/job-%{JOB_ID}-%{TASK_ID}.o.txt" \
        --stderr="${hq_outdir}/job-%{JOB_ID}-%{TASK_ID}.e.txt" \
        -- \
        "${basedir}/hqtask.sh" \
        "${machine}" \
        "${tool}" \
        "${ecf_file}" \
        "${output_path}" \
        "${uniform_clusters_domains}" \
        "${element_type}" \
        "${domains_x}" \
        "${domains_y}" \
        "${domains_z}" \
        "${elements_x}" \
        "${elements_y}" \
        "${elements_z}" \
        "${dual_operator}"
}





if [ "${machine}" == "karolina" ]; then
    sbatch_command_normal="sbatch -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/slurmjob_${machine}.sh"
    sbatch_command_exp="sbatch --partition=qgpu_exp --time=1:00:00 -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/slurmjob_karolina.sh"
    sbatch_command_array="sbatch --array=1-10 -o ${slurm_outdir}/slurm-%A-%a.out -e ${slurm_outdir}/slurm-%A-%a.err ${basedir}/slurmjob_karolina.sh"
elif [ "${machine}" == "lumi" ]; then
    sbatch_command_normal="sbatch -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/slurmjob_${machine}.sh"
    sbatch_command_array="sbatch --array=1-10 -o ${slurm_outdir}/slurm-%A-%a.out -e ${slurm_outdir}/slurm-%A-%a.err ${basedir}/slurmjob_lumi.sh"
elif [ "${machine}" == "e4red" ]; then
    sbatch_command_normal="sbatch -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/slurmjob_${machine}.sh"
elif [[ "${machine}" == "spr"* ]]; then
    sbatch_command_normal="sbatch -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/slurmjob_${machine}.sh"
fi

if [ "${machine}" == "karolina" ]; then
    echo "sbatch command for later, normal job:"
    echo "    ${sbatch_command_normal}"
    echo "sbatch command for later, exp job:"
    echo "    ${sbatch_command_exp}"
    echo "sbatch command for later, array job:"
    echo "    ${sbatch_command_array}"
elif [ "${machine}" == "lumi" ]; then
    echo "sbatch command for later, normal job:"
    echo "    ${sbatch_command_normal}"
    echo "sbatch command for later, array job:"
    echo "    ${sbatch_command_array}"
elif [ "${machine}" == "e4red" ]; then
    echo "sbatch command for later, normal job:"
    echo "    ${sbatch_command_normal}"
elif [ "${machine}" == "tiber" ]; then
    echo "run workers using:"
    echo "    ./benchmarks/dualop_gpu_options/tiber_run_hqworkers.sh"
    echo "inside tmux"
elif [[ "${machine}" == "spr"* ]]; then
    echo "sbatch command for later, normal job:"
    echo "    ${sbatch_command_normal}"
fi




array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf" "${basedir}/espreso.linear_elasticity_2d.ecf")
array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf" "${basedir}/espreso.linear_elasticity_3d.ecf")
array_element_type_2d=(TRIANGLE3 TRIANGLE6)
array_element_type_3d=(TETRA4 TETRA10)
if [ "${machine}" == "karolina" ]; then
    #                   16 32 64 128 256 512 1024 2048
    array_domains_x_2d=( 4  8  8  16  16  32   32  64)
    array_domains_y_2d=( 4  4  8   8  16  16   32  32)
    array_domains_z_2d=( 1  1  1   1   1   1    1   1)
    array_domains_x_3d=( 4  4  4   8   8   8   16  16)
    array_domains_y_3d=( 2  4  4   4   8   8    8  16)
    array_domains_z_3d=( 2  2  4   4   4   8    8   8)
elif [ "${machine}" == "lumi" ]; then
    # approx:           14 28 56 112 224 448 896 1892
    array_domains_x_2d=( 4  7  8  14  16  28  32  56)
    array_domains_y_2d=( 3  4  7   8  14  16  28  32)
    array_domains_z_2d=( 1  1  1   1   1   1   1   1)
    array_domains_x_3d=( 3  4  4   7   8   8  14  16)
    array_domains_y_3d=( 2  3  4   4   7   8   8  14)
    array_domains_z_3d=( 2  2  3   4   4   7   8   8)
elif [ "${machine}" == "e4red" ]; then
    # approx:           24 48 96 192 384 768 1536 3072
    array_domains_x_2d=( 6  8 12  16  24  32   48   64)
    array_domains_y_2d=( 4  6  8  12  16  24   32   48)
    array_domains_z_2d=( 1  1  1   1   1   1   1     1)
    array_domains_x_3d=( 4  4  6   8   8  12  16    16)
    array_domains_y_3d=( 3  4  4   6   8   8  12    16)
    array_domains_z_3d=( 2  3  4   4   6   8   8    12)
elif [ "${machine}" == "tiber" ]; then
    #                   24 48 96 192 384 768 1536 3072
    array_domains_x_2d=( 6  8 12  16  24  32   48  64)
    array_domains_y_2d=( 4  6  8  12  16  24   32  48)
    array_domains_z_2d=( 1  1  1   1   1   1    1   1)
    array_domains_x_3d=( 4  4  6   8   8  12   16  16)
    array_domains_y_3d=( 3  4  4   6   8   8   12  16)
    array_domains_z_3d=( 2  3  4   4   6   8    8  12)
elif [[ "${machine}" == "spr"* ]]; then
    # approx:           12 24 48 96 192 384 768 1536
    array_domains_x_2d=( 4  6  8 12  16  24  32   48)
    array_domains_y_2d=( 3  4  6  8  12  16  24   32)
    array_domains_z_2d=( 1  1  1  1   1   1   1    1)
    array_domains_x_3d=( 3  4  4  6   8   8  12   16)
    array_domains_y_3d=( 2  3  4  4   6   8   8   12)
    array_domains_z_3d=( 2  2  3  4   4   6   8    8)
fi
# indexes:              0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
array_elements_x_2d=(1024 724 512 362 256 181 128  90  64  45  32  23  16  11   8   6)
array_elements_y_2d=(1024 724 512 362 256 181 128  90  64  45  32  23  16  11   8   6)
array_elements_z_2d=(   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1)
array_elements_x_3d=( 101  80  64  50  40  32  25  20  16  13  10   8   6   5   4   3)
array_elements_y_3d=( 101  80  64  50  40  32  25  20  16  13  10   8   6   5   4   3)
array_elements_z_3d=( 101  80  64  50  40  32  25  20  16  13  10   8   6   5   4   3)
array_ecf_file=()
array_element_type=()
array_elements_x=()
array_elements_y=()
array_elements_z=()
array_domains_x=()
array_domains_y=()
array_domains_z=()





uniform_clusters_domains="TRUE"
for dim in 2 3
do
    if [ "${dim}" == "2" ]; then array_ecf_file=("${array_ecf_file_2d[@]}"); else array_ecf_file=("${array_ecf_file_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_element_type=("${array_element_type_2d[@]}"); else array_element_type=("${array_element_type_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_x=("${array_domains_x_2d[@]}"); else array_domains_x=("${array_domains_x_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_y=("${array_domains_y_2d[@]}"); else array_domains_y=("${array_domains_y_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_z=("${array_domains_z_2d[@]}"); else array_domains_z=("${array_domains_z_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_elements_x=("${array_elements_x_2d[@]}"); else array_elements_x=("${array_elements_x_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_elements_y=("${array_elements_y_2d[@]}"); else array_elements_y=("${array_elements_y_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_elements_z=("${array_elements_z_2d[@]}"); else array_elements_z=("${array_elements_z_3d[@]}"); fi
    array_domains_size=${#array_domains_x[@]}
    array_elements_size=${#array_elements_x[@]}

    for ecf_file in "${array_ecf_file[@]}"
    do
        for element_type in "${array_element_type[@]}"
        do
            array_elements_starting_index=20
            if [[ "${ecf_file}" == *"heat_transfer_2d"*     && "${element_type}" == "TRIANGLE3" ]]; then array_elements_starting_index=1; fi
            if [[ "${ecf_file}" == *"heat_transfer_2d"*     && "${element_type}" == "TRIANGLE6" ]]; then array_elements_starting_index=3; fi
            if [[ "${ecf_file}" == *"linear_elasticity_2d"* && "${element_type}" == "TRIANGLE3" ]]; then array_elements_starting_index=2; fi
            if [[ "${ecf_file}" == *"linear_elasticity_2d"* && "${element_type}" == "TRIANGLE6" ]]; then array_elements_starting_index=4; fi
            if [[ "${ecf_file}" == *"heat_transfer_3d"*     && "${element_type}" == "TETRA4"    ]]; then array_elements_starting_index=4; fi
            if [[ "${ecf_file}" == *"heat_transfer_3d"*     && "${element_type}" == "TETRA10"   ]]; then array_elements_starting_index=6; fi
            if [[ "${ecf_file}" == *"linear_elasticity_3d"* && "${element_type}" == "TETRA4"    ]]; then array_elements_starting_index=6; fi
            if [[ "${ecf_file}" == *"linear_elasticity_3d"* && "${element_type}" == "TETRA10"   ]]; then array_elements_starting_index=8; fi

            index=-1
            while true
            do
                index=$((${index}+1))

                array_elements_index=$((${array_elements_starting_index}+${index}))
                array_domains_index=${index}
                if [ "${array_domains_index}" -ge "${array_domains_size}" ]; then array_domains_index=$((${array_domains_size}-1)); fi

                if [ "${array_elements_index}" -ge "${array_elements_size}" ]; then break; fi

                elements_x="${array_elements_x[$array_elements_index]}"
                elements_y="${array_elements_y[$array_elements_index]}"
                elements_z="${array_elements_z[$array_elements_index]}"
                domains_x="${array_domains_x[$array_domains_index]}"
                domains_y="${array_domains_y[$array_domains_index]}"
                domains_z="${array_domains_z[$array_domains_index]}"

                for dual_operator in ${dualoperators}
                do
                    submit_hq_job
                done
            done
        done
    done
done







if [ "${machine}" == "karolina" ]; then
    echo "sbatch command for later, normal job:"
    echo "    ${sbatch_command_normal}"
    echo "sbatch command for later, exp job:"
    echo "    ${sbatch_command_exp}"
    echo "sbatch command for later, array job:"
    echo "    ${sbatch_command_array}"
elif [ "${machine}" == "lumi" ]; then
    echo "sbatch command for later, normal job:"
    echo "    ${sbatch_command_normal}"
    echo "sbatch command for later, array job:"
    echo "    ${sbatch_command_array}"
elif [ "${machine}" == "e4red" ]; then
    echo "sbatch command for later, normal job:"
    echo "    ${sbatch_command_normal}"
elif [ "${machine}" == "tiber" ]; then
    echo "run workers using:"
    echo "    ./benchmarks/dualop_gpu_options/tiber_run_hqworkers.sh"
    echo "inside tmux"
elif [[ "${machine}" == "spr"* ]]; then
    echo "sbatch command for later, normal job:"
    echo "    ${sbatch_command_normal}"
fi

echo "number of jobs: ${id_num}"
