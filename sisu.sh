#!/bin/bash

round()
{
echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
};

module switch PrgEnv-cray/5.2.56 PrgEnv-intel/5.2.56
module switch intel/14.0.4.211 intel/15.0.2.164
module unload cray-libsci
module load gcc/4.9.3
module load cray-tpsl-64

export LC_CTYPE=""

if [ "$#" -ne 1 ]; then
  echo "  Use one of the following commands:"
  echo "    ./sisu.sh configure"
  echo "    ./sisu.sh build"
  echo "    ./sisu.sh run"
  echo "    ./sisu.sh clean"
  echo "    ./sisu.sh distclean"
  echo ""
  echo "  'configure' sets all parameters for the compiler. It is mandatory to run this command at the first time."
  echo "  'build' makes the runable application."
  echo "  'run' starts the computation. The result is put to Paraview input file 'mesh.vtk'."
  echo "  'clean' removes all files created by build process."
  echo "  'distclean' removes all files."
fi

if [ "$1" = "configure" ]; then
  ./waf configure --cray
fi

if [ "$1" = "build" ]; then
  ./waf install -v
fi

if [ "$1" = "clean" ]; then
  ./waf uninstall
  ./waf clean
fi

if [ "$1" = "distclean" ]; then
  ./waf uninstall
  ./waf distclean
fi

#         HEXA8 HEXA20 TETRA4 TETRA10 PRISMA6 PRISMA15 PYRAMID5 PYRAMID13
el_type=(   0     1      2       3       4       5        6         7)

if [ "$1" = "run" ]; then

  export MKL_NUM_THREADS=23
  export OMP_NUM_THREADS=23
  export SOLVER_NUM_THREADS=23
  export PAR_NUM_THREADS=23
  export CILK_NWORKERS=23

  export PARDISOLICMESSAGE=1
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.
  export LC_CTYPE=""

  #               OM OK OK
  #               0   1   2   3   4   5   6   7   8   9
  dom_size=(      8  15  15  15  15  15  15  15  15  15 )
  clustt_size_x=( 8   9   9   9   9   9   9   9   9   9 )
  clustt_size_y=( 8   9   9   9   9   9   9   9   9   9 )
  clustt_size_z=( 9   9   9   9   9   9   9   9   9   9 )

  clusters_x=(    2   3   4   5   6   7   8   9   10   1 )
  clusters_y=(    2   3   4   5   6   7   8   9   10   1 )
  clusters_z=(    2   3   4   5   6   7   8   9   10   1 )

  corners=(       0   0   0   0   0   0   0   0   0   0 )

    #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/composer_xe_2015.2.164/compiler/lib/intel64/:/opt/gcc/4.8.2/snos/lib64/
    #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs

    example_dir="examples/meshgenerator/"
    example="cube_elasticity_fixed_bottom.txt"

    qsub_command_0="#!/bin/bash;"
    qsub_command_0+="export MKL_NUM_THREADS=23;"
    qsub_command_0+="export OMP_NUM_THREADS=23;"
    qsub_command_0+="export SOLVER_NUM_THREADS=23;"
    qsub_command_0+="export PAR_NUM_THREADS=23;"
    qsub_command_0+="export CILK_NWORKERS=23;"
    qsub_command_0+="export PARDISOLICMESSAGE=1;"
    qsub_command_0+="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.;"
    qsub_command_0+="export LC_CTYPE=;"


div=5
K=(  180 120 90 72 60 51 45 40 36 33 30 28 26 24 23 )
cc=(   2   3  4  5  6  7  8  9 10 11 12 13 14 15 16 )
CC=(   1   1  1  1  1  1  1  1  1  1  1  1  1  1  1 )

for div in 4 # 2.5 2.55 2.6 2.65 2.7 2.75 2.8
do

for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do

dt=${K[${i}]}
#d=$(round $dt/$div 0);
dx=$(echo "$dt/$div" | bc -l)
d=`echo "$dt $div" | awk '{printf "%.0f\n", $1/$2}'`
#printf -v $d "%.0f" "$dx"

c=0

x=${cc[${i}]}
y=$x
z=$x

X=${CC[${i}]}
Y=$X
Z=$X



 # for i in 0 # 6 5 0 1 2 3 4
 # do
 #   d=${dom_size[${i}]}
 #   c=${corners[${i}]}

 #   x=${clustt_size_x[${i}]}
 #   y=${clustt_size_y[${i}]}
 #   z=${clustt_size_z[${i}]}

 #   X=${clusters_x[${i}]}
 #   Y=${clusters_y[${i}]}
 #   Z=${clusters_z[${i}]}

    jobname=espreso
    mpiranks=$(( X * Y * Z ))
    account=SERVICE
    actualTime=$( date +%y%m%d_%H:%M:%S )
    log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime.log
    out_dir=ESP-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime
    qsub_command=$qsub_command_0
    qsub_command+="date | tee -a $log_file;"


    $WRKDIR
    echo "Config: dom_size=  $d | cluster_size = $x:$y:$z | clusters = $X:$Y:$Z  "
    mkdir $WRKDIR/$out_dir
    cp -R    ~/espreso/libs/ 	    $WRKDIR/$out_dir
    mkdir -p $WRKDIR/$out_dir/examples
    cp -R    ~/espreso/$example_dir $WRKDIR/$out_dir/$example_dir
    cp       ~/espreso/espreso      $WRKDIR/$out_dir
    cp       ~/espreso/sisu.sh      $WRKDIR/$out_dir

    qsub_command+="cd $WRKDIR/$out_dir;"

    #                                          ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}   | tee -a $log_file
    #mpirun -bind-to-none -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}   | tee -a $log_file
    #mpirun -bind-to none -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}               # | tee -a $log_file


    qsub_command+="aprun -cc none -N 1 -d 23 -n $(( X * Y * Z )) ./espreso $example_dir/$example 0 ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"
    #qsub_command+="aprun -cc none -N 8 -d 3 -n $(( X * Y * Z )) ./espreso $example_dir/$example 0 ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"


    mkdir ~/espreso/logs/$div
    qsub_command+="cp $log_file ~/espreso/logs/$div/"



   echo $qsub_command | tr ";" "\n"
   #echo $qsub_command | tr ";" "\n" | \
   #sbatch -N $(( (X * Y * Z) ))  -p test_large --gid 2000190
   #sbatch -N $(( (X * Y * Z)/8 ))  -p test_large --gid 2000190

  done
 done

fi


if [ "$1" = "debug" ]; then
  module load valgrind/3.9.0-impi
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs
  export OMP_NUM_THREADS=1
  export MKL_PARDISO_OOC_MAX_CORE_SIZE=3500
  export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000

  for i in 0 #1 2 3
  do

    log_file=LOG-1:1:1-2:2:2-5:5:5.log

    date | tee $log_file

    echo "Config: dom_size = 5 | cluster_size = 2:2:2 | clusters = 1:1:1"

    date | tee -a $log_file

    valgrind ./espreso ${el_type[${i}]} 1 1 1 2 2 2 5 5 5 | tee -a $log_file
  done
fi
