#!/bin/bash

machine="karolina"
# machine="mn5"
# machine="lumi"
# machine="tiber"

env_command="source env/it4i.karolina.gcc.cuda.mkl.ss.sh legacy"
# env_command="source env/bsc.mn5.gcc.cuda.mkl.ss.sh legacy"

phase=1





taskid=0

function create_task {

    tasknum=$(printf "%06d" ${taskid})
    taskdir="${tasksdir}/task_${tasknum}"
    mkdir -p "${taskdir}"
    configfile="${taskdir}/config.sh"
    taskid=$((taskid+1))

    (
        echo "#!/bin/bash"
        echo
        echo "${env_command}"
        echo
        echo "export ESPRESO_BENCHMARK_phase=\"${phase}\""
        echo "export ESPRESO_BENCHMARK_machine=\"${machine}\""
        echo "export ESPRESO_BENCHMARK_ecf_file=\"${ecf_file}\""
        echo "export ESPRESO_BENCHMARK_dual_operator=\"${dual_operator}\""
        echo "export ESPRESO_BENCHMARK_element_type=\"${element_type}\""
        echo "export ESPRESO_BENCHMARK_domains_x=\"${domains_x}\""
        echo "export ESPRESO_BENCHMARK_domains_y=\"${domains_y}\""
        echo "export ESPRESO_BENCHMARK_domains_z=\"${domains_z}\""
        echo "export ESPRESO_BENCHMARK_elements_x=\"${elements_x}\""
        echo "export ESPRESO_BENCHMARK_elements_y=\"${elements_y}\""
        echo "export ESPRESO_BENCHMARK_elements_z=\"${elements_z}\""
        echo
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_print_parameters=\"${print_parameters}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_outer_timers=\"${outer_timers}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_inner_timers=\"${inner_timers}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_mainloop_update_split=\"${mainloop_update_split}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_gpu_wait_after_mainloop_update=\"${gpu_wait_after_mainloop_update}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_parallel_set=\"${parallel_set}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_parallel_update=\"${parallel_update}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_parallel_apply=\"${parallel_apply}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_order_F=\"${order_F}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_order_X=\"${order_X}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_order_L=\"${order_L}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_strategy=\"${trsm_strategy}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_partition_algorithm=\"${trsm_partition_algorithm}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_partition_parameter=\"${trsm_partition_parameter}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitrhs_factor_order_sp=\"${trsm_splitrhs_factor_order_sp}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitrhs_factor_order_dn=\"${trsm_splitrhs_factor_order_dn}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitrhs_spdn_criteria=\"${trsm_splitrhs_spdn_criteria}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitrhs_spdn_param=\"${trsm_splitrhs_spdn_param}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_trsm_factor_spdn=\"${trsm_splitfactor_trsm_factor_spdn}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_trsm_factor_order=\"${trsm_splitfactor_trsm_factor_order}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_factor_order_sp=\"${trsm_splitfactor_gemm_factor_order_sp}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_factor_order_dn=\"${trsm_splitfactor_gemm_factor_order_dn}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_factor_prune=\"${trsm_splitfactor_gemm_factor_prune}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_spdn_criteria=\"${trsm_splitfactor_gemm_spdn_criteria}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_trsm_splitfactor_gemm_spdn_param=\"${trsm_splitfactor_gemm_spdn_param}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_herk_strategy=\"${herk_strategy}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_herk_partition_algorithm=\"${herk_partition_algorithm}\""
        echo "export ESPRESO_DUALOPSCTRIA_CONFIG_herk_partition_parameter=\"${herk_partition_parameter}\""
    ) > "${configfile}"
    chmod +x "${configfile}"
}





basedir="benchmarks/dualop_sctria_options"

if [ ! -d "${basedir}" ]
then
    echo "must be run from espreso root directory"
    exit 7
fi

datestr="$(date +%Y%m%d_%H%M%S)"

rundir="${basedir}/runs/${machine}_${datestr}"
mkdir -p "${rundir}"

tasksdir="${rundir}/tasks"
mkdir -p "${tasksdir}"

echo "${machine}" > "${rundir}/machine.txt"

echo "${phase}" > "${rundir}/phase.txt"







# array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf" "${basedir}/espreso.linear_elasticity_2d.ecf")
# array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf" "${basedir}/espreso.linear_elasticity_3d.ecf")
array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf")
array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf")
# array_element_type_2d=(TRIANGLE3 TRIANGLE6)
# array_element_type_3d=(TETRA4 TETRA10)
array_element_type_2d=(TRIANGLE3)
array_element_type_3d=(TETRA4)
if [ "${machine}" == "karolina" ]; then
    #                   16 32 64 128 256 512 1024 2048
    array_domains_x_2d=( 4  8  8  16  16  32   32  64)
    array_domains_y_2d=( 4  4  8   8  16  16   32  32)
    array_domains_z_2d=( 1  1  1   1   1   1    1   1)
    array_domains_x_3d=( 4  4  4   8   8   8   16  16)
    array_domains_y_3d=( 2  4  4   4   8   8    8  16)
    array_domains_z_3d=( 2  2  4   4   4   8    8   8)
elif [ "${machine}" == "mn5" ]; then
    #                   20 40 80 160 320 640 1280 2560
    array_domains_x_2d=( 5  8 10  16  20  32   40  64)
    array_domains_y_2d=( 4  5  8  10  16  20   32  40)
    array_domains_z_2d=( 1  1  1   1   1   1    1   1)
    array_domains_x_3d=( 5  5  5   8   8  10   16  16)
    array_domains_y_3d=( 2  4  4   5   8   8   10  16)
    array_domains_z_3d=( 2  2  4   4   5   8    8  10)
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







print_parameters="1"
outer_timers="1"
inner_timers="0"
mainloop_update_split="S"
gpu_wait_after_mainloop_update="1"
parallel_set="1"
parallel_update="1"
parallel_apply="1"
order_F="R"
order_L="C"

trsm_partition_algorithm="U"
herk_partition_algorithm="U"



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

                if [ "${phase}" == "1" ]
                then
                    #################################################################
                    ### phase 1, select discrete parameters (mainly matrix order) ###
                    #################################################################
                    # numer of tasks: 27x144 = 3888, at most 81 node-hours
                    # fix:
                    #   trsm_partition_parameter
                    #   herk_partition_parameter
                    #   trsm_splitrhs_spdn_param
                    #   trsm_splitfactor_gemm_spdn_param
                    #   herk_strategy (only influence is order_X, and I assume there is none or trsm has greater impact)
                    # for each:
                    #   dual_operator
                    #   trsm_strategy
                    #   trsm_splitrhs_spdn_criteria
                    #   trsm_splitfactor_trsm_factor_spdn
                    #   trsm_splitfactor_gemm_spdn_criteria
                    #   trsm_splitfactor_gemm_factor_prune
                    # find the best combination of:
                    #   order_X
                    #   trsm_splitrhs_factor_order_sp
                    #   trsm_splitrhs_factor_order_dn
                    #   trsm_splitfactor_trsm_factor_order
                    #   trsm_splitfactor_gemm_factor_order_sp
                    #   trsm_splitfactor_gemm_factor_order_dn
                    trsm_partition_parameter="10"
                    herk_partition_parameter="5"
                    herk_strategy="T"
                    trsm_splitrhs_spdn_param="0"
                    trsm_splitfactor_gemm_spdn_param="0"
                    for dual_operator in EXPLICIT_SCTRIA EXPLICIT_SCTRIA_GPU
                    do
                        for order_X in R C
                        do
                            for trsm_strategy in R F
                            do
                                if [ "${trsm_strategy}" == "R" ]
                                then
                                    trsm_splitfactor_trsm_factor_spdn="_"
                                    trsm_splitfactor_gemm_spdn_criteria="_"
                                    trsm_splitfactor_gemm_factor_prune="_"
                                    trsm_splitfactor_trsm_factor_order="_"
                                    trsm_splitfactor_gemm_factor_order_sp="_"
                                    trsm_splitfactor_gemm_factor_order_dn="_"
                                    for trsm_splitrhs_spdn_criteria in S D
                                    do
                                        if [ "${trsm_splitrhs_spdn_criteria}" == "S" ]
                                        then
                                            trsm_splitrhs_factor_order_dn="_"
                                            for trsm_splitrhs_factor_order_sp in R C
                                            do
                                                create_task
                                            done
                                        fi
                                        if [ "${trsm_splitrhs_spdn_criteria}" == "D" ]
                                        then
                                            trsm_splitrhs_factor_order_sp="_"
                                            for trsm_splitrhs_factor_order_dn in R C
                                            do
                                                create_task
                                            done
                                        fi
                                    done
                                fi
                                if [ "${trsm_strategy}" == "F" ]
                                then
                                    trsm_splitrhs_spdn_criteria="_"
                                    trsm_splitrhs_factor_order_sp="_"
                                    trsm_splitrhs_factor_order_dn="_"
                                    for trsm_splitfactor_trsm_factor_spdn in S D
                                    do
                                        for trsm_splitfactor_trsm_factor_order in R C
                                        do
                                            for trsm_splitfactor_gemm_factor_prune in N R
                                            do
                                                for trsm_splitfactor_gemm_spdn_criteria in S D
                                                do
                                                    if [ "${trsm_splitfactor_gemm_spdn_criteria}" == "S" ]
                                                    then
                                                        trsm_splitfactor_gemm_factor_order_dn="_"
                                                        for trsm_splitfactor_gemm_factor_order_sp in R C
                                                        do
                                                            create_task
                                                        done
                                                    fi
                                                    if [ "${trsm_splitfactor_gemm_spdn_criteria}" == "D" ]
                                                    then
                                                        trsm_splitfactor_gemm_factor_order_sp="_"
                                                        for trsm_splitfactor_gemm_factor_order_dn in R C
                                                        do
                                                            create_task
                                                        done
                                                    fi
                                                done
                                            done
                                        done
                                    done
                                fi
                            done
                        done
                    done
                fi


                if [ "${phase}" == "2" ]
                then
                    ########################################################################
                    ### phase 2, select splitfactor_gemm_spdn_param for densiTy criteria ###
                    ########################################################################
                    echo 222
                fi



                if [ "${phase}" == "3" ]
                then
                    ######################################################
                    ### phase 3, partition parameter for trsm and herk ###
                    ######################################################
                    echo 333
                fi



            done
        done
    done
done





echo "submit tasks to HQ using:"
echo "    ./benchmarks/dualop_sctria_options/submit_tasks.sh ${rundir}"
