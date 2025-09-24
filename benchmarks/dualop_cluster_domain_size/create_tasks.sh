#!/bin/bash

# assuming MN5



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
        echo "source env/bsc.mn5.gcc.ompi.mkl.cuda.sh legacy"
        echo
        echo "export ESPRESO_BENCHMARK_ecf_file=\"${ecf_file}\""
        echo
        echo "export ESP_OUTPUT_PATH=${taskdir}"
        echo "export ESP_UNIFORM_DECOMPOSITION=${uniform_decomposition}"
        echo "export ESP_ELEMENT_TYPE=${element_type}"
        echo "export ESP_DOMAINS_X=${domains_x}"
        echo "export ESP_DOMAINS_Y=${domains_y}"
        echo "export ESP_DOMAINS_Z=${domains_z}"
        echo "export ESP_ELEMENTS_X=${elements_x}"
        echo "export ESP_ELEMENTS_Y=${elements_y}"
        echo "export ESP_ELEMENTS_Z=${elements_z}"
        echo "export ESP_FETI_METHOD=${feti_method}"
        echo "export ESP_DUAL_OPERATOR=${dualop}"
        echo "export ESP_PRECONDITIONER=${preconditioner}"
        echo "export ESP_PROJECTOR=${projector}"
    ) > "${configfile}"
    chmod +x "${configfile}"
}





basedir="benchmarks/dualop_cluster_domain_size"

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







array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf" "${basedir}/espreso.linear_elasticity_2d.ecf")
array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf" "${basedir}/espreso.linear_elasticity_3d.ecf")
# array_ecf_file_2d=("${basedir}/espreso.heat_transfer_2d.ecf")
# array_ecf_file_3d=("${basedir}/espreso.heat_transfer_3d.ecf")
array_element_type_2d=(TRIANGLE3 TRIANGLE6 SQUARE4 SQUARE8)
array_element_type_3d=(TETRA4 TETRA10 HEXA8 HEXA20)
# array_element_type_2d=(TRIANGLE3 TRIANGLE6)
# array_element_type_3d=(TETRA4 TETRA10)
# array_element_type_2d=(TRIANGLE3)
# array_element_type_3d=(TETRA4)

# max 40M elements
array_total_cluster_elems_2d_x=(20 20 40 40 80  80 160 160 320 320 640  640 1280 1280 2560 2560 5120 5120)
array_total_cluster_elems_2d_y=(16 32 32 64 64 128 128 256 256 512 512 1024 1024 2048 2048 4096 4096 8192)
array_total_cluster_elems_2d_z=( 1  1  1  1  1   1   1   1   1   1   1    1    1    1    1    1    1    1)
# max 5M elements
array_total_cluster_elems_3d_x=(5 5 10 10 10 20 20 20 40 40 40 80  80  80 160 160)
array_total_cluster_elems_3d_y=(8 8  8 16 16 16 32 32 32 64 64 64 128 128 128 256)
array_total_cluster_elems_3d_z=(4 8  8  8 16 16 16 32 32 32 64 64  64 128 128 128)

# max 5k domains
array_num_domains_2d_x=(1 2 2 5 5 5 10 10 20 20 40 40 80)
array_num_domains_2d_y=(1 1 2 2 4 8  8 16 16 32 32 64 64)
array_num_domains_2d_z=(1 1 1 1 1 1  1  1  1  1  1  1  1)
# max 5k domains
array_num_domains_3d_x=(1 1 1 5 5 5 5 5 5 10 10 10 20)
array_num_domains_3d_y=(1 2 2 2 2 4 4 8 8  8 16 16 16)
array_num_domains_3d_z=(1 1 2 1 2 2 4 4 8  8  8 16 16)

array_feti_methods=(TOTAL_FETI HYBRID_FETI)
array_dual_operators=(IMPLICIT_GENERALSPARSESOLVER_CPU EXPLICIT_GENERALSCHUR_CPU EXPLICIT_GENERALSCHUR_GPU)



uniform_decomposition=true
preconditioner="NONE"
projector="ORTHOGONAL"

for dim in 2 3
do
    if [ "${dim}" == "2" ]; then array_ecf_file=("${array_ecf_file_2d[@]}"); else array_ecf_file=("${array_ecf_file_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_element_type=("${array_element_type_2d[@]}"); else array_element_type=("${array_element_type_3d[@]}"); fi
    if [ "${dim}" == "2" ]; then array_total_cluster_elems_x=("${array_total_cluster_elems_2d_x[@]}"); else array_total_cluster_elems_x=("${array_total_cluster_elems_3d_x[@]}"); fi
    if [ "${dim}" == "2" ]; then array_total_cluster_elems_y=("${array_total_cluster_elems_2d_y[@]}"); else array_total_cluster_elems_y=("${array_total_cluster_elems_3d_y[@]}"); fi
    if [ "${dim}" == "2" ]; then array_total_cluster_elems_z=("${array_total_cluster_elems_2d_z[@]}"); else array_total_cluster_elems_z=("${array_total_cluster_elems_3d_z[@]}"); fi
    if [ "${dim}" == "2" ]; then array_num_domains_x=("${array_num_domains_2d_x[@]}"); else array_num_domains_x=("${array_num_domains_3d_x[@]}"); fi
    if [ "${dim}" == "2" ]; then array_num_domains_y=("${array_num_domains_2d_y[@]}"); else array_num_domains_y=("${array_num_domains_3d_y[@]}"); fi
    if [ "${dim}" == "2" ]; then array_num_domains_z=("${array_num_domains_2d_z[@]}"); else array_num_domains_z=("${array_num_domains_3d_z[@]}"); fi
    array_cluster_elems_size=${#array_total_cluster_elems_x[@]}
    array_num_domains_size=${#array_num_domains_x[@]}

    for ((i=0;i<array_cluster_elems_size;i++))
    do
        for ((j=0;j<array_num_domains_size;j++))
        do
            total_cluster_elems_x="${array_total_cluster_elems_x[$i]}"
            total_cluster_elems_y="${array_total_cluster_elems_y[$i]}"
            total_cluster_elems_z="${array_total_cluster_elems_z[$i]}"
            domains_x="${array_num_domains_x[$j]}"
            domains_y="${array_num_domains_y[$j]}"
            domains_z="${array_num_domains_z[$j]}"
            elements_x=$((total_cluster_elems_x/domains_x))
            elements_y=$((total_cluster_elems_y/domains_y))
            elements_z=$((total_cluster_elems_z/domains_z))

            # min domain size 4x2 or 2x2x2 (so there is some interior node)
            total_elements=$((elements_x*elements_y*elements_z))
            if [ "${total_elements}" -lt 8 ]
            then
                continue
            fi

            if [ $((elements_x*domains_x)) != $total_cluster_elems_x ]; then echo "bad division $total_cluster_elems_x $domains_x"; fi
            if [ $((elements_y*domains_y)) != $total_cluster_elems_y ]; then echo "bad division $total_cluster_elems_y $domains_y"; fi
            if [ $((elements_z*domains_z)) != $total_cluster_elems_z ]; then echo "bad division $total_cluster_elems_z $domains_z"; fi

            for ecf_file in "${array_ecf_file[@]}"
            do
                for element_type in "${array_element_type[@]}"
                do
                    for feti_method in "${array_feti_methods[@]}"
                    do
                        for dualop in "${array_dual_operators[@]}"
                        do
                            create_task
                        done
                    done
                done
            done
        done
    done
done



echo "created ${taskid} tasks"

echo "submit tasks to HQ using:"
echo "    ./benchmarks/dualop_cluster_domain_size/submit_tasks.sh ${rundir}"
