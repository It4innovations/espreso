#!/bin/bash

# assuming Karolina



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
        echo "export OMP_NUM_THREADS=1"
        echo
        echo "export ESPRESO_BENCHMARK_ecf_file=\"${ecf_file}\""
        echo
        echo "export ECF_OUTPUT_PATH=${taskdir}"
        echo "export ECF_ELEMENT_TYPE=${element_type}"
        echo "export ECF_DOMAINS_X=${domains_x}"
        echo "export ECF_DOMAINS_Y=${domains_y}"
        echo "export ECF_DOMAINS_Z=${domains_z}"
        echo "export ECF_ELEMENTS_X=${elements_x}"
        echo "export ECF_ELEMENTS_Y=${elements_y}"
        echo "export ECF_ELEMENTS_Z=${elements_z}"
        echo "export ECF_DUALOPERATOR=${dualop}"
        echo "export ECF_SCHUR_IMPL=${schur}"
        echo "export ECF_SPARSE_SOLVER_IMPL=${spsolver}"
        echo "export ECF_EXPLICIT_APPLY_WHERE=${apply_where}"
    ) > "${configfile}"
    chmod +x "${configfile}"
}





basedir="benchmarks/sparse_solvers_compare"

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







# array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf" "${basedir}/espreso.linear_elasticity_2d.ecf")
# array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf" "${basedir}/espreso.linear_elasticity_3d.ecf")
array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf")
array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf")
# array_element_type_2d=(TRIANGLE3 TRIANGLE6)
# array_element_type_3d=(TETRA4 TETRA10)
array_element_type_2d=(TRIANGLE3)
array_element_type_3d=(TETRA4)

# #                   16 32 64 128 256 512 1024 2048
# array_domains_x_2d=( 4  8  8  16  16  32   32  64)
# array_domains_y_2d=( 4  4  8   8  16  16   32  32)
# array_domains_z_2d=( 1  1  1   1   1   1    1   1)
# array_domains_x_3d=( 4  4  4   8   8   8   16  16)
# array_domains_y_3d=( 2  4  4   4   8   8    8  16)
# array_domains_z_3d=( 2  2  4   4   4   8    8   8)

# decreasing the number of domains 4x because I run singlecore, to save time
#                    4  8 16 32 64 128 256 512
array_domains_x_2d=( 2  4  4  8  8  16  16  32)
array_domains_y_2d=( 2  2  4  4  8   8  16  16)
array_domains_z_2d=( 1  1  1  1  1   1   1   1)
array_domains_x_3d=( 2  2  4  4  4   8   8   8)
array_domains_y_3d=( 2  2  2  4  4   4   8   8)
array_domains_z_3d=( 1  2  2  2  4   4   4   8)

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



dualop_schur_spsolver_triplets=(
    "EXPLICIT_GENERALSCHUR_CPU         AUTO           AUTO         GPU"
    "EXPLICIT_GENERALSCHUR_CPU         MKLPARDISO     AUTO         CPU"
    "EXPLICIT_GENERALSCHUR_CPU         MUMPS          AUTO         CPU"
    "EXPLICIT_GENERALSCHUR_CPU         TRIANGULAR     AUTO         CPU"
    "EXPLICIT_GENERALSCHUR_CPU         SPARSE_SOLVER  SUITESPARSE  CPU"
    "EXPLICIT_GENERALSCHUR_CPU         SPARSE_SOLVER  STRUMPACK    CPU"
    "EXPLICIT_GENERALSCHUR_CPU         SPARSE_SOLVER  PASTIX       CPU"
    "IMPLICIT_GENERALSPARSESOLVER_CPU  AUTO           MKLPARDISO   AUTO"
    "IMPLICIT_GENERALSPARSESOLVER_CPU  AUTO           SUITESPARSE  AUTO"
    "IMPLICIT_GENERALSPARSESOLVER_CPU  AUTO           MUMPS        AUTO"
    "IMPLICIT_GENERALSPARSESOLVER_CPU  AUTO           STRUMPACK    AUTO"
    "IMPLICIT_GENERALSPARSESOLVER_CPU  AUTO           PASTIX       AUTO"
)



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

                for triplet in "${dualop_schur_spsolver_triplets[@]}"
                do
                    set -- $triplet
                    dualop=$1
                    schur=$2
                    spsolver=$3
                    apply_where=$4

                    create_task

                done

            done
        done
    done
done



echo "created ${taskid} tasks"

echo "submit tasks to HQ using:"
echo "    ./benchmarks/sparse_solvers_compare/submit_tasks.sh ${rundir}"
