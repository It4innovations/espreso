#!/bin/bash

source /mnt/proj3/open-18-15/xkadlu01/espreso/scripts/utils.sh

if [ "$1" = "barbora" ]
then
    if [ "$2" = "cpu" ]
    then
        trap 'cleanup_freq_barbora_cn' TERM
        trap 'cleanup_freq_barbora_cn' INT
    else
        trap 'cleanup_freq_barbora_acn' TERM
        trap 'cleanup_freq_barbora_acn' INT
    fi
else
    if [ "$2" = "cpu" ]
    then
        trap 'cleanup_freq_karolina_cn' TERM
        trap 'cleanup_freq_karolina_cn' INT
    else
        trap 'cleanup_freq_karolina_acn' TERM
        trap 'cleanup_freq_karolina_acn' INT
    fi
fi


NODE=$(uname -n | cut -d '.' -f 1)
VERSION_STRING=$(date +"%Y-%m-%d_%Hh%Mm%Ss")

SCR=/scratch/project/open-18-15/$USER/${PBS_JOBNAME}_${VERSION_STRING}
rm -rf SCR
mkdir -p $SCR ; cd $SCR || exit

echo "Start copiyng"
cp -r /mnt/proj3/open-18-15/xkadlu01/espreso/{env,tests,src,include,waf,wscript} ./
echo "End   copiyng"

if [ "$1" = "barbora" ]
then
    if [ "$2" = "cpu" ]
    then
        start_freq_barbora_cn
    else
        start_freq_barbora_acn
    fi
else
    if [ "$2" = "cpu" ]
    then
        start_freq_karolina_cn
    else
        start_freq_karolina_acn
    fi
fi

# res=`find  tests/heatTransfer/ -name espreso.ecf | grep 2D`
# res="${res}
# `find  tests/linearElasticity/ -name espreso.ecf | grep plane`"
# readarray -td$'\n' tasks2D <<<"$res";

# res=`find  tests/heatTransfer/ -name espreso.ecf | grep 3D`
# res="${res}
# `find  tests/linearElasticity/ -name espreso.ecf | grep volume`"
# readarray -td$'\n' tasks3D <<<"$res";

declare -a tasks2D=("tests/heatTransfer/const/conductivity/isotropic2D/espreso.ecf" 
                    "tests/heatTransfer/const/conductivity/anisotropic2D/espreso.ecf")

declare -a tasks3D=("tests/heatTransfer/const/conductivity/isotropic3D/espreso.ecf" 
                    "tests/heatTransfer/const/conductivity/anisotropic3D/espreso.ecf")

declare -a elements3D=("TETRA4" "PYRAMID5" "PRISMA6" "HEXA8" "TETRA10"  "PYRAMID13"  "PRISMA15"  "HEXA20")
declare -a elements_i3D=("8" "8" "8" "8" "8" "8" "8" "8")
declare -a elements_j3D=("8" "8" "8" "8" "8" "8" "8" "8")
declare -a elements_k3D=("8" "8" "8" "8" "8" "8" "8" "8")

declare -a elements2D=("TRIANGLE3" "SQUARE4" "TRIANGLE6" "SQUARE8")
declare -a elements_i2D=("16" "16" "16" "16")
declare -a elements_j2D=("16" "16" "16" "16")

declare -a COMPILERS=("INTEL" "GCC")


# get length of an array
elements3D_len=${#elements3D[@]}
elements2D_len=${#elements2D[@]}

output_file=SC_paper_data_${1}_${2}_${NODE}_${COMPILER}_${VERSION_STRING}.csv
echo "${CSV_HEADER}" > $output_file

for COMPILER in "${COMPILERS[@]}"
do
    if [ "$COMPILER" = "INTEL" ]
    then
        ml purge
        if [ "$1" = "barbora" ]
        then
            source env/modules.barbora.icpc
        else
            source env/modules.karolina.icpc
        fi
        export KMP_AFFINITY=verbose,norespect,granularity=thread,compact
        FLAGS=--cxxflags="-march=native -mtune=native"
    elif [ "$COMPILER" = "GCC" ]
    then
        ml purge
        if [ "$1" = "barbora" ]
        then
            source env/modules.barbora.gcc
        else
            source env/modules.karolina.gcc
        fi

        export OMP_DISPLAY_ENV=TRUE
        export OMP_DISPLAY_AFFINITY=TRUE

        export OMP_PLACES=cores
        export OMP_PROC_BIND=close
        FLAGS=--cxxflags="-march=native -mtune=native"
    else
        echo "Incorrect compiler. Available compilers are INTEL, GCC"
        exit 1
    fi

    ./waf configure "$FLAGS"
    ./waf
    if [ $1 = "karolina" ]
    then
        THREADS=64
    else
        if [ $2 = "cpu" ]
        then
            THREADS=18
        else
            THREADS=12
        fi
    fi
    

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

            ./build/espreso -c $task ${elements3D[$eid]} 1 1 1  1 1 ${THREADS}  $i $j $k > opt3D.tmp
            rm -r ./results

            cat opt3D.tmp | grep -a SCALING: > res.tmp
            echo "$1;$2;$COMPILER;$task;${elements3D[$eid]};${i}x${j}x${k};3D;$THREADS;$(sed  -r  's/.*[^0-9]+\ ([0-9]+\.*[0-9]*)[^0-9]*/\1/'  res.tmp | sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file
            rm opt3D.tmp res.tmp
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

            ./build/espreso -c $task ${elements2D[$eid]} 1 1  1 ${THREADS}  $i $j > opt2D.tmp
            rm -r ./results


            cat opt2D.tmp | grep -a SCALING: > res.tmp
            echo "$1;$2;$COMPILER;$task;${elements2D[$eid]};${i}x${j};2D;$THREADS;$(sed  -r  's/.*[^0-9]+\ ([0-9]+\.*[0-9]*)[^0-9]*/\1/'  res.tmp | sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file
            rm opt2D.tmp res.tmp
        done
    done
done
echo ${VERSION_STRING} >> $output_file

cp $output_file  /mnt/proj3/open-18-15/xkadlu01/espreso/

if [ "$1" = "barbora" ]
then
    if [ "$2" = "cpu" ]
    then
        finish_freq_barbora_cn
    else
        finish_freq_barbora_acn
    fi
else
    if [ "$2" = "cpu" ]
    then
        finish_freq_karolina_cn
    else
        finish_freq_karolina_acn
    fi
fi