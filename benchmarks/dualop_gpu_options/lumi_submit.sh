#!/bin/bash

machine="lumi"
tool="rocm"

basedir="benchmarks/dualop_gpu_options"
if [ ! -f "${basedir}/${machine}_submit.sh" ]
then
    echo "The script must be run from the espreso root folder"
    exit 1
fi

slurm_outdir_root="${basedir}/slurm_outerr"
if [[ ! -d "${slurm_outdir_root}" ]]
then
    echo "Slurm output directory '${slurm_outdir_root}' does not exist"
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
if [[ ! -d "${hq_bin_dir}" ]]
then
    (
        echo "HyperQueue not downloaded, downloading"
        mkdir -p "${hq_bin_dir}"
        cd "${hq_bin_dir}"
        wget https://github.com/It4innovations/hyperqueue/releases/download/v0.19.0/hq-v0.19.0-linux-x64.tar.gz
        tar -xf hq-v0.19.0-linux-x64.tar.gz
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
    output_path="${espreso_outdir}/${benchmarking_for}/${id_str}"
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
        "${concurrency_set}" \
        "${concurrency_update}" \
        "${concurrency_apply}" \
        "${trs1_factor_storage}" \
        "${trs2_factor_storage}" \
        "${trs1_solve_type}" \
        "${trs2_solve_type}" \
        "${trsm_rhs_sol_order}" \
        "${path_if_hermitian}" \
        "${f_sharing_if_hermitian}" \
        "${apply_scatter_gather_where}" \
        "${transpose_where}" \
        "${dualoperator}"
}







sbatch_command_normal="sbatch -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/lumi_slurmjob.sh"
sbatch_command_array="sbatch --array=1-100 -o ${slurm_outdir}/slurm-%A-%a.out -e ${slurm_outdir}/slurm-%A-%a.err ${basedir}/lumi_slurmjob.sh"

echo
echo "sbatch command for later, normal job:"
echo "    ${sbatch_command_normal}"
echo
echo "sbatch command for later, array job:"
echo "    ${sbatch_command_array}"
echo

# ${sbatch_command_array}







array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf" "${basedir}/espreso.linear_elasticity_2d.ecf")
array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf" "${basedir}/espreso.linear_elasticity_3d.ecf")
array_element_type_2d=(TRIANGLE3 TRIANGLE6)
array_element_type_3d=(TETRA4 TETRA10)
array_domains_x_2d=(4 6  8 11 16 23 32 45)
array_domains_y_2d=(4 5  8 11 16 22 32 45)
array_domains_z_2d=(1 1  1  1  1  1  1  1)
array_domains_x_3d=(4 3  4  5  7  8 10 13)
array_domains_y_3d=(2 3  4  5  6  8 10 13)
array_domains_z_3d=(2 3  4  5  6  8 10 12)
array_domains_size=${#array_domains_x_2d[@]}
# indexes:              0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
array_elements_x_2d=(1024 724 512 362 256 181 128  90  64  45  32  23  16  11   8   6)
array_elements_y_2d=(1024 724 512 362 256 181 128  90  64  45  32  23  16  11   8   6)
array_elements_z_2d=(   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1)
array_elements_x_3d=( 101  80  64  50  40  32  25  20  16  13  10   8   6   5   4   3)
array_elements_y_3d=( 101  80  64  50  40  32  25  20  16  13  10   8   6   5   4   3)
array_elements_z_3d=( 101  80  64  50  40  32  25  20  16  13  10   8   6   5   4   3)
array_elements_size=${#array_elements_x_2d[@]}
array_ecf_file=()
array_element_type=()
array_elements_x=()
array_elements_y=()
array_elements_z=()
array_domains_x=()
array_domains_y=()
array_domains_z=()





f_sharing_if_hermitian="SHARED"
transpose_where="CPU"
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

                dualoperator="EXPLICIT_GPU"
                benchmarking_for="setupdate"
                apply_scatter_gather_where="GPU"
                for concurrency in SEQ_WAIT SEQ_CONTINUE PARALLEL
                do
                    concurrency_set="${concurrency}"
                    concurrency_update="${concurrency}"
                    concurrency_apply="${concurrency}"

                    for trs1_factor_storage in SPARSE DENSE
                    do
                        for trsm_rhs_sol_order in ROW_MAJOR COL_MAJOR
                        do
                            for trs1_solve_type in L LHH
                            do
                                path_if_hermitian="HERK"
                                trs2_factor_storage="SPARSE"
                                trs2_solve_type="U"
                                submit_hq_job

                                path_if_hermitian="TRSM"
                                for trs2_factor_storage in SPARSE DENSE
                                do
                                    for trs2_solve_type in U UHH
                                    do
                                        submit_hq_job
                                    done
                                done
                            done
                        done
                    done
                done

                dualoperator="EXPLICIT_GPU"
                benchmarking_for="apply"
                trs1_factor_storage="SPARSE"
                trs1_solve_type="LHH"
                path_if_hermitian="HERK"
                trs2_factor_storage="SPARSE"
                trs2_solve_type="U"
                trsm_rhs_sol_order="ROW_MAJOR"
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

                dualoperator="IMPLICIT_GPU"
                benchmarking_for="setupdateapply"
                apply_scatter_gather_where="GPU"
                path_if_hermitian="TRSM"
                trsm_rhs_sol_order="ROW_MAJOR"
                for concurrency in SEQ_WAIT SEQ_CONTINUE PARALLEL
                do
                    concurrency_set="${concurrency}"
                    concurrency_update="${concurrency}"
                    concurrency_apply="${concurrency}"

                    for trs1_factor_storage in SPARSE DENSE
                    do
                        for trs1_solve_type in L LHH
                        do
                            for trs2_factor_storage in SPARSE DENSE
                            do
                                for trs2_solve_type in U UHH
                                do
                                    submit_hq_job
                                done
                            done
                        done
                    done
                done



                # for dualoperator in EXPLICIT_GPU IMPLICIT_GPU
                # do
                #     benchmarking_for="setupdateapply"
                #     apply_scatter_gather_where="AUTO"
                #     concurrency_set="AUTO"
                #     concurrency_update="AUTO"
                #     concurrency_apply="AUTO"
                #     path_if_hermitian="AUTO"
                #     trs1_factor_storage="AUTO"
                #     trs2_factor_storage="AUTO"
                #     trsm_rhs_sol_order="AUTO"
                #     trs1_solve_type="AUTO"
                #     trs2_solve_type="AUTO"
                #     submit_hq_job
                # done

            done
        done
    done
done





echo
echo "sbatch command for later, normal job:"
echo "    ${sbatch_command_normal}"
echo
echo "sbatch command for later, array job:"
echo "    ${sbatch_command_array}"
echo

echo "number of jobs: ${id_num}"
