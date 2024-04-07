#!/bin/bash

machine="lumi"
tool="rocm"

basedir="benchmarks/dualop_cpuimpl_vs_gpuexpl"
if [ ! -f "${basedir}/${machine}_submit.sh" ]
then
    echo "The script must be run from the espreso root folder"
    exit 1
fi

slurm_outdir_root="${basedir}/slurm_outerr"
if [[ ! -d "${slurm_outdir_root}"   ||   ! ( -L "${slurm_outdir_root}" && -d "${slurm_outdir_root}" ) ]]
then
    echo "Slurm output directory '${slurm_outdir_root}' does not exist"
    echo "Please create it"
    exit 2
fi

espreso_outdir_root="${basedir}/espreso_results"
if [[ ! -d "${espreso_outdir_root}"   ||   ! ( -L "${espreso_outdir_root}" && -d "${espreso_outdir_root}" ) ]]
then
    echo "Espreso output directory '${espreso_outdir_root}' does not exist"
    echo "Please create it"
    exit 3
fi

hq_outdir_root="${basedir}/hq_outerr"
if [[ ! -d "${hq_outdir_root}"   ||   ! ( -L "${hq_outdir_root}" && -d "${hq_outdir_root}" ) ]]
then
    echo "HyperQueue output directory '${hq_outdir_root}' does not exist"
    echo "Please create it"
    exit 4
fi

hq_bin_dir="dependencies/HyperQueue"
if [[ ! -d "${hq_bin_dir}" ]]
then
    (
        echo "HyperQueue not downloaded, downloading"
        mkdir -p "${hq_bin_dir}"
        cd "${hq_bin_dir}"
        wget https://github.com/It4innovations/hyperqueue/releases/download/v0.18.0/hq-v0.18.0-linux-x64.tar.gz
        tar -xf hq-v0.18.0-linux-x64.tar.gz
    )
fi

HQBIN="${hq_bin_dir}/hq"
if ! "${HQBIN}" server info > /dev/null 2> /dev/null
then
    echo "HyperQueue server is not running"
    echo "Start the server on this login node using '${HQBIN} server start' inside tmux"
    exit 5
fi



datestr="$(date +%Y%m%d_%H%M%S)"

hq_outdir="${hq_outdir_root}/${machine}_${tool}_${datestr}"
mkdir -p "${hq_outdir}"

espreso_outdir="${espreso_outdir_root}/${machine}_${tool}_${datestr}"
mkdir -p "${espreso_outdir}"

slurm_outdir="${slurm_outdir_root}/${machine}_${tool}_${datestr}"
mkdir -p "${slurm_outdir}"






id_num=0
function submit_hq_job
{
    id_num=$((${id_num}+1))
    id_str=$(printf '%06d' ${id_num})
    output_path="${espreso_outdir}/${id_str}"
    echo "Submitting task ${id_str}"
    "${HQBIN}" submit \
        --time-request=6m \
        --cpus='7 compact' \
        --stdout="${hq_outdir}/job-%{JOB_ID}-%{TASK_ID}.o.txt" \
        --stderr="${hq_outdir}/job-%{JOB_ID}-%{TASK_ID}.e.txt" \
        -- \
        "${basedir}/lumi_hqtask.sh" \
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







sbatch_command_normal="sbatch -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/lumi_slurmjob.sh"
sbatch_command_array="sbatch --array=1-10 -o ${slurm_outdir}/slurm-%A-%a.out -e ${slurm_outdir}/slurm-%A-%a.err ${basedir}/lumi_slurmjob.sh"

echo "sbatch command for later, normal job:"
echo "    ${sbatch_command_normal}"
echo "sbatch command for later, array job:"
echo "    ${sbatch_command_array}"

# ${sbatch_command_array}







array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf" "${basedir}/espreso.linear_elasticity_2d.ecf")
array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf" "${basedir}/espreso.linear_elasticity_3d.ecf")
array_element_type_2d=(TRIANGLE3 TRIANGLE6)
array_element_type_3d=(TETRA4 TETRA10)
array_domains_x_2d=(32 23 16 11 8)
array_domains_y_2d=(32 22 16 11 8)
array_domains_z_2d=( 1  1  1  1 1)
array_domains_x_3d=(10 8 7 5 4)
array_domains_y_3d=(10 8 6 5 4)
array_domains_z_3d=(10 8 6 5 4)
array_ecf_file=()
array_element_type=()
array_elements_x=()
array_elements_y=()
array_elements_z=()
array_domains_x=()
array_domains_y=()
array_domains_z=()





for dim in 2 3
do
    if [ "${dim}" == "2" ]; then array_ecf_file=("${array_ecf_file_2d[@]}"); else array_ecf_file=("${array_ecf_file_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_element_type=("${array_element_type_2d[@]}"); else array_element_type=("${array_element_type_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_x=("${array_domains_x_2d[@]}"); else array_domains_x=("${array_domains_x_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_y=("${array_domains_y_2d[@]}"); else array_domains_y=("${array_domains_y_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_z=("${array_domains_z_2d[@]}"); else array_domains_z=("${array_domains_z_3d[@]}"); fi

    for ecf_file in "${array_ecf_file[@]}"
    do
        for element_type in "${array_element_type[@]}"
        do
            if [[ "${ecf_file}" == *"heat_transfer_2d"* ]]; then
                array_elements_x=(8 16 32 64 128)
                array_elements_y=(8 16 32 64 128)
                array_elements_z=(1  1  1  1   1)
            fi
            if [[ "${ecf_file}" == *"linear_elasticity_2d"* ]]; then
                array_elements_x=(8 16 32 64 128)
                array_elements_y=(8 16 32 64 128)
                array_elements_z=(1  1  1  1   1)
            fi
            if [[ "${ecf_file}" == *"heat_transfer_3d"* && "${element_type}" == "TETRA4" ]]; then
                array_elements_x=(4 6 10 16 25)
                array_elements_y=(4 6 10 16 25)
                array_elements_z=(4 6 10 16 25)
            fi
            if [[ "${ecf_file}" == *"heat_transfer_3d"* && "${element_type}" == "TETRA10" ]]; then
                array_elements_x=(3 4 6 10 16)
                array_elements_y=(3 4 6 10 16)
                array_elements_z=(3 4 6 10 16)
            fi
            if [[ "${ecf_file}" == *"linear_elasticity_3d"* && "${element_type}" == "TETRA4" ]]; then
                array_elements_x=(3 4 6 10 16)
                array_elements_y=(3 4 6 10 16)
                array_elements_z=(3 4 6 10 16)
            fi
            if [[ "${ecf_file}" == *"linear_elasticity_3d"* && "${element_type}" == "TETRA10" ]]; then
                array_elements_x=(3 4 5 6 8)
                array_elements_y=(3 4 5 6 8)
                array_elements_z=(3 4 5 6 8)
            fi

            for i in ${!array_domains_x[@]}
            do
                elements_x="${array_elements_x[$i]}"
                elements_y="${array_elements_y[$i]}"
                elements_z="${array_elements_z[$i]}"
                domains_x="${array_domains_x[$i]}"
                domains_y="${array_domains_y[$i]}"
                domains_z="${array_domains_z[$i]}"

                for uniform_clusters_domains in TRUE FALSE
                do
                    for dual_operator in IMPLICIT EXPLICIT_GPU
                    do
                        submit_hq_job
                    done
                done
            done
        done
    done
done





echo "sbatch command for later, normal job:"
echo "    ${sbatch_command_normal}"
echo "sbatch command for later, array job:"
echo "    ${sbatch_command_array}"
