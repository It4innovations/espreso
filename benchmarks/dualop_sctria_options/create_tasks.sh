#!/bin/bash

machine="karolina"
# machine="mn5"
# machine="lumi"
# machine="tiber"

env_command="source env/it4i.karolina.gcc.cuda.mkl.ss.sh legacy"
# env_command="source env/bsc.mn5.gcc.cuda.mkl.ss.sh legacy"



if [ "$#" -lt 1 ]
then
    echo "not enough arguments"
    exit 2
fi
phase="${1}"





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

rundir="${basedir}/runs/${machine}_phase${phase}_${datestr}"
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
                    # numer of tasks: 27x144 = 3888 = at most 324 gpu-hours
                    # fix:
                    #   trsm_partition_algorithm
                    #   herk_partition_algorithm
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
                    trsm_partition_algorithm="U"
                    herk_partition_algorithm="U"
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
                    ###########################################
                    ### phase 2, select trsm spdn and prune ###
                    ###########################################
                    # numer of tasks: 27x80 = 2160 = at most 180 gpu-hours
                    # auto-select what we have already determined:
                    #   order_X
                    #   trsm_splitrhs_factor_order_sp
                    #   trsm_splitrhs_factor_order_dn
                    #   trsm_splitfactor_trsm_factor_order
                    #   trsm_splitfactor_gemm_factor_order_sp
                    #   trsm_splitfactor_gemm_factor_order_dn
                    # fix:
                    #   trsm_partition_algorithm
                    #   herk_partition_algorithm
                    #   trsm_splitrhs_spdn_param
                    #   trsm_splitfactor_gemm_spdn_param
                    #   herk_strategy
                    #   herk_partition_parameter
                    # for each:
                    #   dual_operator
                    #   trsm_strategy
                    #   trsm_partition_parameter
                    # select the best combination of:
                    #   trsm_splitrhs_spdn_criteria
                    #   trsm_splitfactor_trsm_factor_spdn
                    #   trsm_splitfactor_gemm_spdn_criteria
                    #   trsm_splitfactor_gemm_factor_prune
                    order_X="_"
                    trsm_splitrhs_factor_order_sp="_"
                    trsm_splitrhs_factor_order_dn="_"
                    trsm_splitfactor_trsm_factor_order="_"
                    trsm_splitfactor_gemm_factor_order_sp="_"
                    trsm_splitfactor_gemm_factor_order_dn="_"
                    trsm_partition_algorithm="U"
                    herk_partition_algorithm="U"
                    trsm_splitrhs_spdn_param="0"
                    trsm_splitfactor_gemm_spdn_param="0"
                    herk_strategy="T"
                    herk_partition_parameter="5"
                    partition_parameters=()
                    partition_parameters_2d_domain=(-20000 -2000 10 100)
                    partition_parameters_3d_domain=(-2000 -200 10 100)
                    partition_parameters_2d_interface=(-500 -50 10 100)
                    partition_parameters_3d_interface=(-1000 -100 10 100)
                    for dual_operator in EXPLICIT_SCTRIA EXPLICIT_SCTRIA_GPU
                    do
                        for trsm_strategy in R F
                        do
                            if [ "${dim}" == "2" ] && [ "${trsm_strategy}" == "R" ]; then partition_parameters=("${partition_parameters_2d_interface[@]}"); fi
                            if [ "${dim}" == "3" ] && [ "${trsm_strategy}" == "R" ]; then partition_parameters=("${partition_parameters_3d_interface[@]}"); fi
                            if [ "${dim}" == "2" ] && [ "${trsm_strategy}" == "F" ]; then partition_parameters=("${partition_parameters_2d_domain[@]}"); fi
                            if [ "${dim}" == "3" ] && [ "${trsm_strategy}" == "F" ]; then partition_parameters=("${partition_parameters_3d_domain[@]}"); fi
                            for trsm_partition_parameter in "${partition_parameters[@]}"
                            do
                                if [ "${trsm_strategy}" == "R" ]
                                then
                                    trsm_splitfactor_trsm_factor_spdn="_"
                                    trsm_splitfactor_gemm_spdn_criteria="_"
                                    trsm_splitfactor_gemm_factor_prune="_"
                                    for trsm_splitrhs_spdn_criteria in S D
                                    do
                                        create_task
                                    done
                                fi
                                if [ "${trsm_strategy}" == "F" ]
                                then
                                    trsm_splitrhs_spdn_criteria="_"
                                    for trsm_splitfactor_trsm_factor_spdn in S D
                                    do
                                        for trsm_splitfactor_gemm_spdn_criteria in S D
                                        do
                                            for trsm_splitfactor_gemm_factor_prune in N R
                                            do
                                                create_task
                                            done
                                        done 
                                    done
                                fi
                            done
                        done
                    done
                fi



                if [ "${phase}" == "3" ]
                then
                    #########################################
                    ### phase 3, trsm partition parameter ###
                    #########################################
                    # numer of tasks: 2244 = at most 187 gpu-hours
                    # auto-select what we have already determined:
                    #   order_X
                    #   trsm_splitrhs_factor_order_sp
                    #   trsm_splitrhs_factor_order_dn
                    #   trsm_splitfactor_trsm_factor_order
                    #   trsm_splitfactor_gemm_factor_order_sp
                    #   trsm_splitfactor_gemm_factor_order_dn
                    #   trsm_splitrhs_spdn_criteria
                    #   trsm_splitfactor_trsm_factor_spdn
                    #   trsm_splitfactor_gemm_spdn_criteria
                    #   trsm_splitfactor_gemm_factor_prune
                    # fix:
                    #   trsm_partition_algorithm
                    #   herk_partition_algorithm
                    #   trsm_splitrhs_spdn_param
                    #   trsm_splitfactor_gemm_spdn_param
                    #   herk_strategy
                    #   herk_partition_parameter
                    # for each:
                    #   dual_operator
                    #   trsm_strategy
                    # select the best:
                    #   trsm_partition_parameter
                    order_X="_"
                    trsm_splitrhs_factor_order_sp="_"
                    trsm_splitrhs_factor_order_dn="_"
                    trsm_splitfactor_trsm_factor_order="_"
                    trsm_splitfactor_gemm_factor_order_sp="_"
                    trsm_splitfactor_gemm_factor_order_dn="_"
                    trsm_splitrhs_spdn_criteria="_"
                    trsm_splitfactor_trsm_factor_spdn="_"
                    trsm_splitfactor_gemm_spdn_criteria="_"
                    trsm_splitfactor_gemm_factor_prune="_"
                    trsm_partition_algorithm="U"
                    herk_partition_algorithm="U"
                    trsm_splitrhs_spdn_param="0"
                    trsm_splitfactor_gemm_spdn_param="0"
                    herk_strategy="T"
                    herk_partition_parameter="5"
                    partition_parameters=()
                    partition_parameters_2d_splitfactor=(-100000 -50000 -20000 -10000 -5000 -2000 -1000 -500 -200 -100 -50 -20 -10 1 2 5 10 20 50 100 200 500 1000)
                    partition_parameters_3d_splitfactor=(-50000 -20000 -10000 -5000 -2000 -1000 -500 -200 -100 -50 -20 -10 1 2 5 10 20 50 100 200 500 1000)
                    partition_parameters_2d_splitrhs=(-5000 -2000 -1000 -500 -200 -100 -50 -20 -10 1 2 5 10 20 50 100 200 500 1000)
                    partition_parameters_3d_splitrhs=(-5000 -2000 -1000 -500 -200 -100 -50 -20 -10 1 2 5 10 20 50 100 200 500 1000)
                    for dual_operator in EXPLICIT_SCTRIA EXPLICIT_SCTRIA_GPU
                    do
                        for trsm_strategy in R F
                        do
                            if [ "${dim}" == "2" ] && [ "${trsm_strategy}" == "R" ]; then partition_parameters=("${partition_parameters_2d_splitrhs[@]}"); fi
                            if [ "${dim}" == "3" ] && [ "${trsm_strategy}" == "R" ]; then partition_parameters=("${partition_parameters_3d_splitrhs[@]}"); fi
                            if [ "${dim}" == "2" ] && [ "${trsm_strategy}" == "F" ]; then partition_parameters=("${partition_parameters_2d_splitfactor[@]}"); fi
                            if [ "${dim}" == "3" ] && [ "${trsm_strategy}" == "F" ]; then partition_parameters=("${partition_parameters_3d_splitfactor[@]}"); fi
                            for trsm_partition_parameter in "${partition_parameters[@]}"
                            do
                                create_task
                            done
                        done
                    done
                fi



                if [ "${phase}" == "4" ]
                then
                    #########################################
                    ### phase 4, herk partition parameter ###
                    #########################################
                    # numer of tasks: 2244 = at most 187 gpu-hours
                    # auto-select what we have already determined:
                    #   order_X
                    #   trsm_splitrhs_factor_order_sp
                    #   trsm_splitrhs_factor_order_dn
                    #   trsm_splitfactor_trsm_factor_order
                    #   trsm_splitfactor_gemm_factor_order_sp
                    #   trsm_splitfactor_gemm_factor_order_dn
                    #   trsm_splitrhs_spdn_criteria
                    #   trsm_splitfactor_trsm_factor_spdn
                    #   trsm_splitfactor_gemm_spdn_criteria
                    #   trsm_splitfactor_gemm_factor_prune
                    #   trsm_partition_parameter
                    # fix:
                    #   trsm_partition_algorithm
                    #   herk_partition_algorithm
                    #   trsm_splitrhs_spdn_param
                    #   trsm_splitfactor_gemm_spdn_param
                    #   trsm_strategy
                    # for each:
                    #   dual_operator
                    #   herk_strategy
                    # select the best:
                    #   herk_partition_parameter
                    order_X="_"
                    trsm_splitrhs_factor_order_sp="_"
                    trsm_splitrhs_factor_order_dn="_"
                    trsm_splitfactor_trsm_factor_order="_"
                    trsm_splitfactor_gemm_factor_order_sp="_"
                    trsm_splitfactor_gemm_factor_order_dn="_"
                    trsm_splitrhs_spdn_criteria="_"
                    trsm_splitfactor_trsm_factor_spdn="_"
                    trsm_splitfactor_gemm_spdn_criteria="_"
                    trsm_splitfactor_gemm_factor_prune="_"
                    trsm_partition_parameter="0"
                    trsm_partition_algorithm="U"
                    herk_partition_algorithm="U"
                    trsm_splitrhs_spdn_param="0"
                    trsm_splitfactor_gemm_spdn_param="0"
                    trsm_strategy="F"
                    partition_parameters=()
                    partition_parameters_2d_squares=(-100000 -50000 -20000 -10000 -5000 -2000 -1000 -500 -200 -100 -50 -20 -10 1 2 5 10 20 50 100 200 500 1000)
                    partition_parameters_3d_squares=(-50000 -20000 -10000 -5000 -2000 -1000 -500 -200 -100 -50 -20 -10 1 2 5 10 20 50 100 200 500 1000)
                    partition_parameters_2d_stairs=(-5000 -2000 -1000 -500 -200 -100 -50 -20 -10 1 2 5 10 20 50 100 200 500 1000)
                    partition_parameters_3d_stairs=(-5000 -2000 -1000 -500 -200 -100 -50 -20 -10 1 2 5 10 20 50 100 200 500 1000)
                    for dual_operator in EXPLICIT_SCTRIA EXPLICIT_SCTRIA_GPU
                    do
                        for herk_strategy in T Q
                        do
                            if [ "${dim}" == "2" ] && [ "${herk_strategy}" == "T" ]; then partition_parameters=("${partition_parameters_2d_stairs[@]}"); fi
                            if [ "${dim}" == "3" ] && [ "${herk_strategy}" == "T" ]; then partition_parameters=("${partition_parameters_3d_stairs[@]}"); fi
                            if [ "${dim}" == "2" ] && [ "${herk_strategy}" == "Q" ]; then partition_parameters=("${partition_parameters_2d_squares[@]}"); fi
                            if [ "${dim}" == "3" ] && [ "${herk_strategy}" == "Q" ]; then partition_parameters=("${partition_parameters_3d_squares[@]}"); fi
                            for herk_partition_parameter in "${partition_parameters[@]}"
                            do
                                create_task
                            done
                        done
                    done
                fi



                if [ "${phase}" == "5" ]
                then
                    #######################################
                    ### phase 5, trsm and herk strategy ###
                    #######################################
                    # numer of tasks: 27x24 = 648 = at most 54 gpu-hours
                    # auto-select what we have already determined:
                    #   order_X
                    #   trsm_splitrhs_factor_order_sp
                    #   trsm_splitrhs_factor_order_dn
                    #   trsm_splitfactor_trsm_factor_order
                    #   trsm_splitfactor_gemm_factor_order_sp
                    #   trsm_splitfactor_gemm_factor_order_dn
                    #   trsm_splitrhs_spdn_criteria
                    #   trsm_splitfactor_trsm_factor_spdn
                    #   trsm_splitfactor_gemm_spdn_criteria
                    #   trsm_partition_parameter
                    #   herk_partition_parameter
                    # fix:
                    #   trsm_splitrhs_spdn_param
                    #   trsm_splitfactor_gemm_spdn_param
                    # for each:
                    #   dual_operator
                    #   trsm_splitfactor_gemm_factor_prune
                    # select the best:
                    #   trsm_partition_algorithm
                    #   herk_partition_algorithm
                    #   herk_strategy
                    #   trsm_strategy
                    order_X="_"
                    trsm_splitrhs_factor_order_sp="_"
                    trsm_splitrhs_factor_order_dn="_"
                    trsm_splitfactor_trsm_factor_order="_"
                    trsm_splitfactor_gemm_factor_order_sp="_"
                    trsm_splitfactor_gemm_factor_order_dn="_"
                    trsm_splitrhs_spdn_criteria="_"
                    trsm_splitfactor_trsm_factor_spdn="_"
                    trsm_splitfactor_gemm_spdn_criteria="_"
                    trsm_partition_parameter="0"
                    herk_partition_parameter="0"
                    trsm_splitrhs_spdn_param="0"
                    trsm_splitfactor_gemm_spdn_param="0"
                    for dual_operator in EXPLICIT_SCTRIA EXPLICIT_SCTRIA_GPU
                    do
                        herk_partition_algorithm="_"
                        herk_strategy="_"
                        for trsm_partition_algorithm in U M
                        do
                            for trsm_strategy in R F
                            do
                                for trsm_splitfactor_gemm_factor_prune in N R
                                do
                                    create_task
                                done
                            done
                        done

                        trsm_partition_algorithm="_"
                        trsm_strategy="_"
                        trsm_splitfactor_gemm_factor_prune="_"
                        for herk_partition_algorithm in U M
                        do
                            for herk_strategy in T Q
                            do
                                create_task
                            done
                        done
                    done
                fi



                if [ "${phase}" == "6" ]
                then
                    ###################################################################
                    ### phase 6, compare with original algorithm (set only 1 chunk) ###
                    ###################################################################
                    # numer of tasks: 27x8 = 216 = at most 18 gpu-hours
                    # auto-select what we have already determined (almost everything):
                    #   order_X
                    #   trsm_splitrhs_factor_order_sp
                    #   trsm_splitrhs_factor_order_dn
                    #   trsm_splitfactor_trsm_factor_order
                    #   trsm_splitfactor_gemm_factor_order_sp
                    #   trsm_splitfactor_gemm_factor_order_dn
                    #   trsm_splitrhs_spdn_criteria
                    #   trsm_splitfactor_trsm_factor_spdn
                    #   trsm_splitfactor_gemm_spdn_criteria
                    #   trsm_splitfactor_gemm_factor_prune
                    #   herk_strategy
                    #   trsm_strategy
                    #   trsm_splitrhs_spdn_param
                    #   trsm_splitfactor_gemm_spdn_param
                    #   trsm_partition_algorithm
                    #   herk_partition_algorithm
                    # for each:
                    #   dual_operator
                    #   mainloop_update_split
                    # observe the difference between:
                    #   trsm_partition_parameter=1 and autoselect
                    #   herk_partition_parameter=1 and autoselect
                    order_X="_"
                    trsm_splitrhs_factor_order_sp="_"
                    trsm_splitrhs_factor_order_dn="_"
                    trsm_splitfactor_trsm_factor_order="_"
                    trsm_splitfactor_gemm_factor_order_sp="_"
                    trsm_splitfactor_gemm_factor_order_dn="_"
                    trsm_splitfactor_gemm_spdn_criteria="_"
                    trsm_splitfactor_gemm_factor_prune="_"
                    herk_strategy="_"
                    trsm_strategy="_"
                    trsm_splitrhs_spdn_param="0"
                    trsm_splitfactor_gemm_spdn_param="0"
                    trsm_partition_algorithm="U"
                    herk_partition_algorithm="U"
                    total_elements=$((${elements_x} * ${elements_y} * ${elements_z}))
                    for dual_operator in EXPLICIT_SCTRIA EXPLICIT_SCTRIA_GPU
                    do
                        for mainloop_update_split in S C
                        do
                            for pp in 0 1
                            do
                                if [ "${pp}" == 1 ] && [ "${dual_operator}" == "EXPLICIT_SCTRIA_GPU" ]
                                then
                                    trsm_splitrhs_spdn_criteria="S"
                                    trsm_splitfactor_trsm_factor_spdn="S"
                                    if [ "${dim}" == "3" ] && [ "${total_elements}" -lt "12000" ]
                                    then
                                        trsm_splitrhs_spdn_criteria="D"
                                        trsm_splitfactor_trsm_factor_spdn="D"
                                    fi
                                else
                                    trsm_splitrhs_spdn_criteria="_"
                                    trsm_splitfactor_trsm_factor_spdn="_"
                                fi

                                trsm_partition_parameter="${pp}"
                                herk_partition_parameter="${pp}"

                                create_task
                            done
                        done
                    done
                fi



                if [ "${phase}" == "7" ]
                then
                    #########################################################################
                    ### phase 7, compare pure trsm and herk kernel times with numchunks=1 ###
                    #########################################################################
                    # numer of tasks: 27x4 = 108 = at most 9 gpu-hours
                    # auto-select what we have already determined (almost everything):
                    #   order_X
                    #   trsm_splitrhs_factor_order_sp
                    #   trsm_splitrhs_factor_order_dn
                    #   trsm_splitfactor_trsm_factor_order
                    #   trsm_splitfactor_gemm_factor_order_sp
                    #   trsm_splitfactor_gemm_factor_order_dn
                    #   trsm_splitrhs_spdn_criteria
                    #   trsm_splitfactor_trsm_factor_spdn
                    #   trsm_splitfactor_gemm_spdn_criteria
                    #   trsm_splitfactor_gemm_factor_prune
                    #   herk_strategy
                    #   trsm_strategy
                    #   trsm_splitrhs_spdn_param
                    #   trsm_splitfactor_gemm_spdn_param
                    #   trsm_partition_algorithm
                    #   herk_partition_algorithm
                    # fix because I need to measure some things:
                    #   parallel_set
                    #   parallel_update
                    #   parallel_apply
                    # for each:
                    #   dual_operator
                    # observe the difference between:
                    #   trsm_partition_parameter=1 and autoselect
                    #   herk_partition_parameter=1 and autoselect
                    order_X="_"
                    trsm_splitrhs_factor_order_sp="_"
                    trsm_splitrhs_factor_order_dn="_"
                    trsm_splitfactor_trsm_factor_order="_"
                    trsm_splitfactor_gemm_factor_order_sp="_"
                    trsm_splitfactor_gemm_factor_order_dn="_"
                    trsm_splitfactor_gemm_spdn_criteria="_"
                    trsm_splitfactor_gemm_factor_prune="_"
                    herk_strategy="_"
                    trsm_strategy="_"
                    trsm_splitrhs_spdn_param="0"
                    trsm_splitfactor_gemm_spdn_param="0"
                    trsm_partition_algorithm="U"
                    herk_partition_algorithm="U"
                    parallel_set="0"
                    parallel_update="0"
                    parallel_apply="0"
                    total_elements=$((${elements_x} * ${elements_y} * ${elements_z}))
                    for dual_operator in EXPLICIT_SCTRIA EXPLICIT_SCTRIA_GPU
                    do
                        for pp in 0 1
                        do
                            if [ "${pp}" == 1 ] && [ "${dual_operator}" == "EXPLICIT_SCTRIA_GPU" ]
                            then
                                trsm_splitrhs_spdn_criteria="S"
                                trsm_splitfactor_trsm_factor_spdn="S"
                                if [ "${dim}" == "3" ] && [ "${total_elements}" -lt "12000" ]
                                then
                                    trsm_splitrhs_spdn_criteria="D"
                                    trsm_splitfactor_trsm_factor_spdn="D"
                                fi
                            else
                                trsm_splitrhs_spdn_criteria="_"
                                trsm_splitfactor_trsm_factor_spdn="_"
                            fi

                            trsm_partition_parameter="${pp}"
                            herk_partition_parameter="${pp}"

                            create_task
                        done
                    done
                fi



                if [ "${phase}" == "8" ]
                then
                    ######################################################################
                    ### phase 8, compare pruning and non-pruning and sparse/dense gemm ###
                    ######################################################################
                    # numer of tasks: 27x8 = 216 = atmost 18 gpu-hours
                    # auto-select what we have already determined (almost everything):
                    #   order_X
                    #   trsm_splitrhs_factor_order_sp
                    #   trsm_splitrhs_factor_order_dn
                    #   trsm_splitfactor_trsm_factor_order
                    #   trsm_splitfactor_gemm_factor_order_sp
                    #   trsm_splitfactor_gemm_factor_order_dn
                    #   trsm_splitrhs_spdn_criteria
                    #   trsm_splitfactor_trsm_factor_spdn
                    #   herk_strategy
                    #   trsm_splitrhs_spdn_param
                    #   trsm_splitfactor_gemm_spdn_param
                    #   trsm_partition_algorithm
                    #   herk_partition_algorithm
                    #   trsm_partition_parameter
                    #   herk_partition_parameter
                    # fix:
                    #   trsm_strategy
                    # for each:
                    #   dual_operator
                    # observe the difference between:
                    #   trsm_splitfactor_gemm_factor_prune
                    #   trsm_splitfactor_gemm_spdn_criteria
                    order_X="_"
                    trsm_splitrhs_factor_order_sp="_"
                    trsm_splitrhs_factor_order_dn="_"
                    trsm_splitfactor_trsm_factor_order="_"
                    trsm_splitfactor_gemm_factor_order_sp="_"
                    trsm_splitfactor_gemm_factor_order_dn="_"
                    trsm_splitrhs_spdn_criteria="_"
                    trsm_splitfactor_trsm_factor_spdn="_"
                    herk_strategy="_"
                    trsm_splitrhs_spdn_param="0"
                    trsm_splitfactor_gemm_spdn_param="0"
                    trsm_partition_algorithm="U"
                    herk_partition_algorithm="U"
                    trsm_partition_parameter="0"
                    herk_partition_parameter="0"
                    trsm_strategy="F"
                    for dual_operator in EXPLICIT_SCTRIA EXPLICIT_SCTRIA_GPU
                    do
                        for trsm_splitfactor_gemm_spdn_criteria in S D
                        do
                            for trsm_splitfactor_gemm_factor_prune in N R
                            do
                                create_task
                            done
                        done
                    done
                fi



            done
        done
    done
done



echo "created ${taskid} tasks"

echo "submit tasks to HQ using:"
echo "    ./benchmarks/dualop_sctria_options/submit_tasks.sh ${rundir}"
