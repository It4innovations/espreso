#!/bin/bash

machine="karolina"

# have to make sure the same version is compiled
tool="cudalegacy"
# tool="cudamodern"

basedir="benchmarks/dualop_gpu_options"
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

module load HyperQueue/0.18.0
if ! hq server info > /dev/null 2> /dev/null
then
    echo "HyperQueue server is not running"
    echo "Start the server on this login node using 'hq server start' inside tmux"
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
    output_path="${espreso_outdir}/${benchmarking_for}/${id_str}"
    echo "Submitting task ${id_str}"
    hq submit \
        --time-request=6m \
        --cpus='16 compact' \
        --stdout="${hq_outdir}/job-%{JOB_ID}-%{TASK_ID}.o.txt" \
        --stderr="${hq_outdir}/job-%{JOB_ID}-%{TASK_ID}.e.txt" \
        -- \
        "${basedir}/karolina_hqtask.sh" \
        "${machine}" \
        "${tool}" \
        "${ecf_file}" \
        "${factor_symmetry_fake}" \
        "${output_path}" \
        "${uniform_clusters_domains}" \
        "${element_type}" \
        "${domains_x}" \
        "${domains_y}" \
        "${domains_z}" \
        "${elements_x}" \
        "${elements_y}" \
        "${elements_z}" \
        "${concurrency_set}" \
        "${concurrency_update}" \
        "${concurrency_apply}" \
        "${trsm1_factor_storage}" \
        "${trsm2_factor_storage}" \
        "${trsm1_solve_type}" \
        "${trsm2_solve_type}" \
        "${trsm_rhs_sol_order}" \
        "${path_if_hermitian}" \
        "${f_sharing_if_hermitian}" \
        "${apply_scatter_gather_where}" \
        "${transpose_where}"
}







sbatch_command_normal="sbatch -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/karolina_slurmjob.sh"
sbatch_command_exp="sbatch --partition=qgpu_exp --time=1:00:00 -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/karolina_slurmjob.sh"
sbatch_command_array="sbatch --array=1-100 -o ${slurm_outdir}/slurm-%A-%a.out -e ${slurm_outdir}/slurm-%A-%a.err ${basedir}/karolina_slurmjob.sh"

${sbatch_command_array}







array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf" "${basedir}/espreso.linear_elasticity_2d.ecf")
array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf" "${basedir}/espreso.linear_elasticity_3d.ecf")
array_element_type_2d=(TRIANGLE3 TRIANGLE6)
array_element_type_3d=(TETRA4 TETRA10)
array_elements_x_2d=(8 16 32 64 128)
array_elements_y_2d=(8 16 32 64 128)
array_elements_z_2d=(1  1  1  1   1)
array_elements_x_3d=(4 6 10 16 25)
array_elements_y_3d=(4 6 10 16 25)
array_elements_z_3d=(4 6 10 16 25)
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





# set and update
benchmarking_for="setupdate"
f_sharing_if_hermitian="SHARED"
apply_scatter_gather_where="GPU"
transpose_where="GPU"
for dim in 2 3
do
    if [ "${dim}" == "2" ]; then array_ecf_file=("${array_ecf_file_2d[@]}"); else array_ecf_file=("${array_ecf_file_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_element_type=("${array_element_type_2d[@]}"); else array_element_type=("${array_element_type_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_elements_x=("${array_elements_x_2d[@]}"); else array_elements_x=("${array_elements_x_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_elements_y=("${array_elements_y_2d[@]}"); else array_elements_y=("${array_elements_y_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_elements_z=("${array_elements_z_2d[@]}"); else array_elements_z=("${array_elements_z_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_x=("${array_domains_x_2d[@]}"); else array_domains_x=("${array_domains_x_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_y=("${array_domains_y_2d[@]}"); else array_domains_y=("${array_domains_y_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_z=("${array_domains_z_2d[@]}"); else array_domains_z=("${array_domains_z_3d[@]}"); fi

    for ecf_file in "${array_ecf_file[@]}"
    do
        for element_type in "${array_element_type[@]}"
        do
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
                    for concurrency in SEQ_WAIT SEQ_CONTINUE PARALLEL
                    do
                        concurrency_set="${concurrency}"
                        concurrency_update="${concurrency}"
                        concurrency_apply="${concurrency}"

                        for factor_symmetry_fake in L U B
                        do
                            for trsm1_factor_storage in SPARSE DENSE
                            do
                                for trsm1_solve_type in L LHH
                                do
                                    for trsm_rhs_sol_order in ROW_MAJOR COL_MAJOR
                                    do
                                        if [[ "${factor_symmetry_fake}" != "B" ]]
                                        then
                                            path_if_hermitian="HERK"
                                            trsm2_factor_storage="SPARSE"
                                            trsm2_solve_type="U"
                                            submit_hq_job
                                        fi

                                        path_if_hermitian="TRSM"
                                        for trsm2_factor_storage in SPARSE DENSE
                                        do
                                            for trsm2_solve_type in U UHH
                                            do
                                                submit_hq_job
                                            done
                                        done
                                    done
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done





# apply
benchmarking_for="apply"
factor_symmetry_fake="U"
trsm1_factor_storage="SPARSE"
trsm1_solve_type="L"
path_if_hermitian="HERK"
trsm2_factor_storage="DENSE"
trsm2_solve_type="U"
trsm_rhs_sol_order="ROW_MAJOR"
f_sharing_if_hermitian="SHARED"
transpose_where="GPU"
for dim in 2 3
do
    if [ "${dim}" == "2" ]; then array_ecf_file=("${array_ecf_file_2d[@]}"); else array_ecf_file=("${array_ecf_file_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_element_type=("${array_element_type_2d[@]}"); else array_element_type=("${array_element_type_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_elements_x=("${array_elements_x_2d[@]}"); else array_elements_x=("${array_elements_x_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_elements_y=("${array_elements_y_2d[@]}"); else array_elements_y=("${array_elements_y_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_elements_z=("${array_elements_z_2d[@]}"); else array_elements_z=("${array_elements_z_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_x=("${array_domains_x_2d[@]}"); else array_domains_x=("${array_domains_x_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_y=("${array_domains_y_2d[@]}"); else array_domains_y=("${array_domains_y_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_domains_z=("${array_domains_z_2d[@]}"); else array_domains_z=("${array_domains_z_3d[@]}"); fi

    for ecf_file in "${array_ecf_file[@]}"
    do
        for element_type in "${array_element_type[@]}"
        do
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
                    for concurrency in SEQ_WAIT SEQ_CONTINUE PARALLEL
                    do
                        concurrency_set="PARALLEL"
                        concurrency_update="PARALLEL"
                        concurrency_apply="${concurrency}"

                        for apply_scatter_gather_where in CPU GPU
                        do
                            submit_hq_job
                        done
                    done
                done
            done
        done
    done
done





echo "sbatch command for later, normal job:"
echo "    ${sbatch_command_normal}"
echo "sbatch command for later, exp job:"
echo "    ${sbatch_command_exp}"
echo "sbatch command for later, array job:"
echo "    ${sbatch_command_array}"
