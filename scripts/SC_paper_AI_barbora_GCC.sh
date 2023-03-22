#!/bin/bash
#PBS -q qcpu
#PBS -N espreso_ai_barbora_GCC
#PBS -l select=1:ncpus=36,walltime=24:00:00,msr=1.4.0
#PBS -A OPEN-22-7

CLUSTER=barbora
COMPILER=GCC

a=$(date +%s)

source /mnt/proj3/open-18-15/xkadlu01/espreso/scripts/utils.sh

SDE_HOME=/home/xkadlu01/Installs/sde-external-9.7.0-2022-05-09-lin/

NODE=$(uname -n | cut -d '.' -f 1)
VERSION_STRING=$(date +"%Y-%m-%d_%Hh%Mm%Ss")

SCR=/scratch/project/open-18-15/$USER/${PBS_JOBNAME}_${VERSION_STRING}
rm -rf SCR
mkdir -p $SCR ; cd $SCR || exit

echo "Start copiyng"
cp -r /mnt/proj3/open-18-15/xkadlu01/espreso/{env,tests,prebuild,SDE-testbed} ./
echo "End   copiyng"

# res=`find  tests/heatTransfer/ -name espreso.ecf | grep 2D`
# res="${res}
# `find  tests/linearElasticity/ -name espreso.ecf | grep plane`"
# readarray -td$'\n' tasks2D <<<"$res";

# res=`find  tests/heatTransfer/ -name espreso.ecf | grep 3D`
# res="${res}
# `find  tests/linearElasticity/ -name espreso.ecf | grep volume`"
# readarray -td$'\n' tasks3D <<<"$res";


declare -a tasks2D=("tests/heatTransfer/const/conductivity/isotropic2D/espreso.ecf"
                    "tests/heatTransfer/const/conductivity/symmetric2D/espreso.ecf"
                    "tests/heatTransfer/const/coordinateSystem/cartesian2D/espreso.ecf"
                    "tests/heatTransfer/const/conductivity/anisotropic2D/espreso.ecf"
                    "tests/heatTransfer/const/coordinateSystem/cylindrical2D/espreso.ecf"
                    "tests/heatTransfer/const/elements/heatSource2D/espreso.ecf"
                    "tests/heatTransfer/const/advection/conductivity/isotropic2D/espreso.ecf"
                    "tests/heatTransfer/const/advection/conductivity/diagonal2D/espreso.ecf"
                    "tests/linearElasticity/const/planeAxisymmetric/elements/acceleration/espreso.ecf"
                    "tests/linearElasticity/const/planeStrain/elements/acceleration/espreso.ecf"
                    "tests/linearElasticity/const/planeStress/elements/acceleration/espreso.ecf"
                    "tests/linearElasticity/const/planeStressWithThickness/elements/acceleration/espreso.ecf"
                    "tests/linearElasticity/const/planeAxisymmetric/elements/angularVelocity/espreso.ecf"
                    "tests/linearElasticity/const/planeStrain/elements/angularVelocity/espreso.ecf"
                    "tests/linearElasticity/const/planeStress/elements/angularVelocity/espreso.ecf"
                    "tests/linearElasticity/const/planeStressWithThickness/elements/angularVelocity/espreso.ecf")

# declare -a tasks3D=("tests/heatTransfer/const/conductivity/isotropic3D/espreso.ecf" 
#                     "tests/heatTransfer/const/conductivity/anisotropic3D/espreso.ecf")

declare -a tasks3D=("tests/heatTransfer/const/conductivity/isotropic3D/espreso.ecf"
                    "tests/heatTransfer/const/conductivity/symmetric3D/espreso.ecf"
                    "tests/heatTransfer/const/coordinateSystem/cartesian3D/espreso.ecf"
                    "tests/heatTransfer/const/conductivity/anisotropic3D/espreso.ecf"
                    "tests/heatTransfer/const/coordinateSystem/cylindrical3D/espreso.ecf"
                    "tests/heatTransfer/const/coordinateSystem/spherical3D/espreso.ecf"
                    "tests/heatTransfer/const/elements/heatSource3D/espreso.ecf"
                    "tests/heatTransfer/const/advection/conductivity/isotropic3D/espreso.ecf"
                    "tests/heatTransfer/const/advection/conductivity/diagonal3D/espreso.ecf"
                    "tests/linearElasticity/const/volume/element/acceleration/espreso.ecf"
                    "tests/linearElasticity/const/volume/element/angularVelocityZ/espreso.ecf"
                    "tests/linearElasticity/const/volume/elasticity/orthotropic/espreso.ecf"
                    "tests/linearElasticity/const/volume/coordinateSystem/cartesian/espreso.ecf"
                    "tests/linearElasticity/const/volume/coordinateSystem/cylindrical/espreso.ecf"
                    "tests/linearElasticity/const/volume/output/stress/espreso.ecf")

declare -a elements3D=("TETRA4" "PYRAMID5" "PRISMA6" "HEXA8" "TETRA10"  "PYRAMID13"  "PRISMA15"  "HEXA20")
declare -a elements_i3D=("8" "8" "8" "8" "8" "8" "8" "8")
declare -a elements_j3D=("8" "8" "8" "8" "8" "8" "8" "8")
declare -a elements_k3D=("8" "8" "8" "8" "8" "8" "8" "8")

declare -a elements2D=("TRIANGLE3" "SQUARE4" "TRIANGLE6" "SQUARE8")
declare -a elements_i2D=("16" "16" "16" "16")
declare -a elements_j2D=("16" "16" "16" "16")

# declare -a COMPILERS=("INTEL" "GCC")
# declare -a CLUSTERS=("barbora" "karolina")
declare -a TYPES=("INHERITANCE" "CONDITIONS")


# get length of an array
elements3D_len=${#elements3D[@]}
elements2D_len=${#elements2D[@]}

output_file=SC_paper_data_AI_${NODE}_${VERSION_STRING}.csv
echo "cluster;compiler;task;element;decomposition;dim;threads;impl;action;WRITE AI; READ AI; OVERALL AI; SUM WRITE AI; SUM READ AI; SUM OVERALL AI; UMSKD SP FLOP; MSKD SP FLOP; UMSKD DP FLOP; MSKD DP FLOP; SUM SP FLOP; SUM DP FLOP; FMA; SUM FMA" > $output_file

# for CLUSTER in "${CLUSTERS[@]}"
# do
#     for COMPILER in "${COMPILERS[@]}"
#     do
        if [ "$COMPILER" = "INTEL" ]
        then
            ml purge
            if [ "$CLUSTER" = "barbora" ]
            then
                source env/modules.barbora.icpc
            else
                source env/modules.karolina.icpc
            fi
        elif [ "$COMPILER" = "GCC" ]
        then
            ml purge
            if [ "$CLUSTER" = "barbora" ]
            then
                source env/modules.barbora.gcc
            else
                source env/modules.karolina.gcc
            fi
        else
            echo "Incorrect compiler. Available compilers are INTEL, GCC"
            exit 1
        fi


        THREADS=1
        

        . env/threading.default ${THREADS}

        for task in "${tasks3D[@]}"
        do
            for (( eid=0; eid<${elements3D_len}; eid++ ));
            do
                i=${elements_i3D[$eid]}
                j=${elements_j3D[$eid]}
                k=${elements_k3D[$eid]}

                echo "Testing $task Decomposition into $i x $j x $k elements"
                # use for loop to read all values and indexes
                for TYPE in "${TYPES[@]}"
                do

                    cd $SCR/SDE-testbed
                    ${SDE_HOME}/sde64 -iform -mix -dyn_mask_profile -start_ssc_mark FACE:repeat -stop_ssc_mark DEAD:repeat -- ../prebuild/${CLUSTER}/${COMPILER}/espreso -c ../${task} ${elements3D[$eid]} 1 1 1  1 1 ${THREADS}  $i $j $k --DRYRUN=1 --LOOP=${TYPE}  >/dev/null
                    ./intel_sde_flops.py 1>../AI.tmp
                    rm -r ./results
                    rm sde-dyn*
                    rm sde-mix*
                    cd $SCR

                    echo "$CLUSTER;$COMPILER;$task;${elements3D[$eid]};${i}x${j}x${k};3D;$THREADS;${TYPE};ASSEMBLE;$(sed -n 's/.*\((approx\.): \)\(.*\)\((EXPERIMENTAL)\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FLOPs: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FMA instructions executed: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file

                    cd $SCR/SDE-testbed
                    ${SDE_HOME}/sde64 -iform -mix -dyn_mask_profile -start_ssc_mark CAFE:repeat -stop_ssc_mark FADE:repeat -- ../prebuild/${CLUSTER}/${COMPILER}/espreso -c ../${task} ${elements3D[$eid]} 1 1 1  1 1 ${THREADS}  $i $j $k --DRYRUN=1 --LOOP=${TYPE}  >/dev/null
                    ./intel_sde_flops.py 1>../AI.tmp
                    rm -r ./results
                    rm sde-dyn*
                    rm sde-mix*
                    cd $SCR

                    echo "$CLUSTER;$COMPILER;$task;${elements3D[$eid]};${i}x${j}x${k};3D;$THREADS;${TYPE};RE-ASSEMBLE;$(sed -n 's/.*\((approx\.): \)\(.*\)\((EXPERIMENTAL)\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FLOPs: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FMA instructions executed: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file

                    cd $SCR/SDE-testbed
                    ${SDE_HOME}/sde64 -iform -mix -dyn_mask_profile -start_ssc_mark FEED:repeat -stop_ssc_mark BEED:repeat -- ../prebuild/${CLUSTER}/${COMPILER}/espreso -c ../${task} ${elements3D[$eid]} 1 1 1  1 1 ${THREADS}  $i $j $k --DRYRUN=1 --LOOP=${TYPE}  >/dev/null
                    ./intel_sde_flops.py 1>../AI.tmp
                    rm -r ./results
                    rm sde-dyn*
                    rm sde-mix*
                    cd $SCR

                    echo "$CLUSTER;$COMPILER;$task;${elements3D[$eid]};${i}x${j}x${k};3D;$THREADS;${TYPE};SOLUTION;$(sed -n 's/.*\((approx\.): \)\(.*\)\((EXPERIMENTAL)\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FLOPs: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FMA instructions executed: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file
                done
            done
        done

        for task in "${tasks2D[@]}"
        do
            for (( eid=0; eid<${elements2D_len}; eid++ ));
            do

                i=${elements_i2D[$eid]}
                j=${elements_j2D[$eid]}

                echo "Testing $task Decomposition into 1 x $i x $j elements"
                # use for loop to read all values and indexes
                for TYPE in "${TYPES[@]}"
                do
                    cd $SCR/SDE-testbed
                    ${SDE_HOME}/sde64 -iform -mix -dyn_mask_profile -start_ssc_mark FACE:repeat -stop_ssc_mark DEAD:repeat -- ../prebuild/${CLUSTER}/${COMPILER}/espreso -c ../${task} ${elements2D[$eid]} 1 1  1 ${THREADS}  $i $j --DRYRUN=1 --LOOP=${TYPE} >/dev/null
                    ./intel_sde_flops.py 1>../AI.tmp
                    rm -r ./results
                    rm sde-dyn*
                    rm sde-mix*
                    cd $SCR

                    echo "$CLUSTER;$COMPILER;$task;${elements2D[$eid]};${i}x${j};2D;$THREADS;${TYPE};ASSEMBLE;$(sed -n 's/.*\((approx\.): \)\(.*\)\((EXPERIMENTAL)\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FLOPs: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FMA instructions executed: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file

                    cd $SCR/SDE-testbed
                    ${SDE_HOME}/sde64 -iform -mix -dyn_mask_profile -start_ssc_mark CAFE:repeat -stop_ssc_mark FADE:repeat -- ../prebuild/${CLUSTER}/${COMPILER}/espreso -c ../${task} ${elements2D[$eid]} 1 1  1 ${THREADS}  $i $j --DRYRUN=1 --LOOP=${TYPE} >/dev/null
                    ./intel_sde_flops.py 1>../AI.tmp
                    rm -r ./results
                    rm sde-dyn*
                    rm sde-mix*
                    cd $SCR

                    echo "$CLUSTER;$COMPILER;$task;${elements2D[$eid]};${i}x${j};2D;$THREADS;${TYPE};RE-ASSEMBLE;$(sed -n 's/.*\((approx\.): \)\(.*\)\((EXPERIMENTAL)\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FLOPs: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FMA instructions executed: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file

                    cd $SCR/SDE-testbed
                    ${SDE_HOME}/sde64 -iform -mix -dyn_mask_profile -start_ssc_mark FEED:repeat -stop_ssc_mark BEED:repeat -- ../prebuild/${CLUSTER}/${COMPILER}/espreso -c ../${task} ${elements2D[$eid]} 1 1  1 ${THREADS}  $i $j --DRYRUN=1 --LOOP=${TYPE} >/dev/null
                    ./intel_sde_flops.py 1>../AI.tmp
                    rm -r ./results
                    rm sde-dyn*
                    rm sde-mix*
                    cd $SCR

                    echo "$CLUSTER;$COMPILER;$task;${elements2D[$eid]};${i}x${j};2D;$THREADS;${TYPE};SOLUTION;$(sed -n 's/.*\((approx\.): \)\(.*\)\((EXPERIMENTAL)\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FLOPs: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/');$(sed -n 's/.*\(FMA instructions executed: \)\(.*\)/\2/p' AI.tmp |  sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file
                
                done
            done
        done
#     done
# done

b=$(date +%s)
diff=$((b-a))
diff=$((diff-3600))

DURATION=$(date +"%H:%M:%S" --date=@${diff})
echo ${VERSION_STRING} >> $output_file
echo ${DURATION} >> $output_file

cp $output_file  /mnt/proj3/open-18-15/xkadlu01/espreso/