#!/bin/bash

CSV_HEADER="run;compiler;element;task;decomposition;threads;implementation;input elements;all elements;thread elements;REPETITIONS;MIN INTERVALS;MAX INTERVALS;PARALLEL REGION;THREAD AVG;THREAD MIN;THREAD MAX;THREAD CP;KERNEL SUM AVG;KERNEL SUM MIN;KERNEL SUM MAX;KERNEL AVG;KERNEL MIN;KERNEL MAX;TOTAL DP OPS;THREAD DP OPS;KERNEL DP OPS;OVERALL GFLOPS; THREAD CP GFLOPS; AVG KERNEL GFLOPS; TOTAL DRAM TFC; AVG THREAD DRAM TFC; KERNEL DRAM TFC; TOTAL DRAM BW; AVG THREAD DRAM BW; KERNEL DRAM BW;"




function cleanup_freq_barbora_cn {

	echo "==== WARNING SCRIPT HAS FINISHED ABNORMALY! ====="
	echo "========= REVERTING TIMING TO DEFAULT ==========="

	source /mnt/proj3/open-18-15/meric_barbora/set_env

	$MERIC_ROOT/tools/energyMeasureStop -e RAPL

	exit 1
}

function start_freq_barbora_cn {
    echo "MERIC BARBORA CPU"
    source /mnt/proj3/open-18-15/meric_barbora/set_env 
    $MERIC_ROOT/tools/energyMeasureStart -c 1600000000 -e RAPL #-u 2400000000 -e RAPL
    sleep 5 
}

function finish_freq_barbora_cn {

	echo "========= SCRIPT HAS FINISHED NORMALY ==========="
    echo "========= REVERTING TIMING TO DEFAULT ==========="

	source /mnt/proj3/open-18-15/meric_barbora/set_env

	$MERIC_ROOT/tools/energyMeasureStop -e RAPL
}

function cleanup_freq_barbora_acn {

	echo "==== WARNING SCRIPT HAS FINISHED ABNORMALY! ====="
	echo "========= REVERTING TIMING TO DEFAULT ==========="

	source /mnt/proj3/open-18-15/meric_barbora/set_env

	$MERIC_ROOT/tools/energyMeasureStop -e RAPL

	exit 1
}

function start_freq_barbora_acn {
    echo "MERIC BARBORA GPU"
    source /mnt/proj3/open-18-15/meric_barbora/set_env 
    $MERIC_ROOT/tools/energyMeasureStart -c 1700000000 -e RAPL #-u 2400000000 -e RAPL
    sleep 5 
}

function finish_freq_barbora_acn {

	echo "========= SCRIPT HAS FINISHED NORMALY ==========="
    echo "========= REVERTING TIMING TO DEFAULT ==========="

	source /mnt/proj3/open-18-15/meric_barbora/set_env

	$MERIC_ROOT/tools/energyMeasureStop -e RAPL
}

function cleanup_freq_karolina_cn {

	echo "==== WARNING SCRIPT HAS FINISHED ABNORMALY! ====="
	echo "========= REVERTING TIMING TO DEFAULT ==========="

	source /mnt/proj3/open-18-15/meric_karolina/set_env_cn

	$MERIC_ROOT/tools/energyMeasureStop

	exit 1
}

function start_freq_karolina_cn {
    echo "MERIC KAROLINA CPU"
    source /mnt/proj3/open-18-15/meric_karolina/set_env_cn

    $MERIC_ROOT/tools/energyMeasureStart -c 2600000000 # -u 2400000000
    sleep 5 
}

function finish_freq_karolina_cn {

	echo "========= SCRIPT HAS FINISHED NORMALY ==========="
    echo "========= REVERTING TIMING TO DEFAULT ==========="

	source /mnt/proj3/open-18-15/meric_karolina/set_env_cn

	$MERIC_ROOT/tools/energyMeasureStop
}

function cleanup_freq_karolina_acn {

	echo "==== WARNING SCRIPT HAS FINISHED ABNORMALY! ====="
	echo "========= REVERTING TIMING TO DEFAULT ==========="

	source /mnt/proj3/open-18-15/meric_karolina/set_env_acn

	$MERIC_ROOT/tools/energyMeasureStop

	exit 1
}

function start_freq_karolina_acn {
    echo "MERIC KAROLINA GPU"
    source /mnt/proj3/open-18-15/meric_karolina/set_env_acn 
    $MERIC_ROOT/tools/energyMeasureStart -c 2600000000 # -u 2400000000
    sleep 5 
}

function finish_freq_karolina_acn {

	echo "========= SCRIPT HAS FINISHED NORMALY ==========="
    echo "========= REVERTING TIMING TO DEFAULT ==========="

	source /mnt/proj3/open-18-15/meric_karolina/set_env_acn

	$MERIC_ROOT/tools/energyMeasureStop
}

function run3D {
    . env/threading.default ${THREADS}
    for (( eid=0; eid<${3Delements_len}; eid++ ));
    do

        i=3Delements_i[$eid]
        j=3Delements_j[$eid]
        k=3Delements_k[$eid]

        echo "Testing $3. Decomposition into 1 x $i x $j elements with size $size"
        # use for loop to read all values and indexes
    
        ./build/espreso -c Ecfs/3D_${task_ids[$t]}.ecf 3Delements[$eid] 1 1 1  1 1 16  $i $j $k $impl $measurement > opt3D.tmp
        rm -r ./results


        cat opt3D.tmp | grep -a SCALING: > res.tmp
        echo "$runnumber;$3;$4;${tasks[$t]};${i_print}x${j_print};$thread;$impl;$size;$relsize;$ptsize;$(sed  -r  's/.*[^0-9]+\ ([0-9]+\.*[0-9]*)[^0-9]*/\1/'  res.tmp | sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file
        rm opt3D.tmp res.tmp
    done
}


function run2D {
    . env/threading.default ${THREADS}
    for (( eid=0; eid<${2Delements_len}; eid++ ));
    do

        i=2Delements_i[$eid]
        j=2Delements_j[$eid]

        echo "Testing $3. Decomposition into 1 x $i x $j elements with size $size"
        # use for loop to read all values and indexes
    
        ./build/espreso -c Ecfs/3D_${task_ids[$t]}.ecf 2Delements[$eid] 1 1 1  1 1 16  $i $j $k $impl $measurement > opt2D.tmp
        rm -r ./results


        cat opt2D.tmp | grep -a SCALING: > res.tmp
        echo "$runnumber;$3;$4;${tasks[$t]};${i_print}x${j_print};$thread;$impl;$size;$relsize;$ptsize;$(sed  -r  's/.*[^0-9]+\ ([0-9]+\.*[0-9]*)[^0-9]*/\1/'  res.tmp | sed  -z 's/\n/;/g;s/;$/\n/')" >> $output_file
        rm opt2D.tmp res.tmp
    done
}