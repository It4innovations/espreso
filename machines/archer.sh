#!/bin/bash

WORKDIR=~/work/ondra
ESPRESODIR=~/ondra/espreso
EXAMPLEDIR=examples/meshgenerator
EXAMPLE=cube_elasticity_fixed_bottom.txt

module switch PrgEnv-cray/5.2.56 PrgEnv-intel/5.2.56
module switch intel/14.0.4.211 intel/15.0.2.164 
module unload cray-libsci
module load gcc/4.9.3
module load cray-tpsl
module list

export LC_CTYPE=""

if [ "$#" -ne 1 ]; then
  echo "  Use one of the following commands:"
  echo "    ./archer.sh configure"
  echo "    ./archer.sh build"
  echo "    ./archer.sh run"
  echo ""
  echo "  'configure' sets all parameters for the compiler. It is mandatory to run this command at the first time."
  echo "  'build' makes the runable application."
  echo "  'run' starts the computation. The result is put to Paraview input file 'mesh.vtk'."
fi

if [ "$1" = "configure" ]; then
  ./waf configure --cray
fi

if [ "$1" = "build" ]; then
  ./waf install -v -j4
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
  #               0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  dom_size=(      15  15  15  15  15  15  15  15  15  15  15  15  15  15  15  15  15  15  15  15)
  clustt_size_x=( 9   9   9   9   9   9   9   9   9   9    9   9   9   9   9   9   9   9   9   9)
  clustt_size_y=( 9   9   9   9   9   9   9   9   9   9    9   9   9   9   9   9   9   9   9   9)
  clustt_size_z=( 9   9   9   9   9   9   9   9   9   9    9   9   9   9   9   9   9   9   9   9)

  clusters_x=(    1   2   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21)
  clusters_y=(    1   2   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21)
  clusters_z=(    1   2   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21)

  corners=(       0   0   0   0   0   0   0   0   0    0   0   0   0   0   0   0   0   0   0   0)
    
    qsub_command_0="#!/bin/bash;"
    qsub_command_0+="export MKL_NUM_THREADS=23;"
    qsub_command_0+="export OMP_NUM_THREADS=23;"
    qsub_command_0+="export SOLVER_NUM_THREADS=23;"
    qsub_command_0+="export PAR_NUM_THREADS=23;"
    qsub_command_0+="export CILK_NWORKERS=23;"
    qsub_command_0+="export PARDISOLICMESSAGE=1;"
    qsub_command_0+="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.;"
    qsub_command_0+="export LC_CTYPE=;"


  for i in 0 # 11 12 # 6 5 0 1 2 3 4
  do
    d=${dom_size[${i}]}
    c=${corners[${i}]}

    x=${clustt_size_x[${i}]}
    y=${clustt_size_y[${i}]}
    z=${clustt_size_z[${i}]}

    X=${clusters_x[${i}]}
    Y=${clusters_y[${i}]}
    Z=${clusters_z[${i}]}

    jobname=espreso
    mpiranks=$(( X * Y * Z ))
    account=pr1u0905
    actualTime=$( date +%y%m%d_%H:%M:%S )
    log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime.log
    out_dir=ESP-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime
    qsub_command=$qsub_command_0
    qsub_command+="date | tee -a $log_file;"

    echo "Config: dom_size=  $d | cluster_size = $x:$y:$z | clusters = $X:$Y:$Z  "
    mkdir    $WORKDIR/$out_dir
    mkdir -p $WORKDIR/$out_dir/$EXAMPLEDIR
    cp -R    $ESPRESODIR/libs/ 	        $WORKDIR/$out_dir
    cp -R    $ESPRESODIR/$EXAMPLEDIR/*  $WORKDIR/$out_dir/$EXAMPLEDIR
    cp       $ESPRESODIR/espreso        $WORKDIR/$out_dir
    cp       $ESPRESODIR/archer.sh      $WORKDIR/$out_dir
    
    qsub_command+="cd $WORKDIR/$out_dir;"
    qsub_command+="ls;"

    #                                          ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}   | tee -a $log_file
    #mpirun -bind-to-none -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}   | tee -a $log_file
    #mpirun -bind-to none -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}               # | tee -a $log_file


    qsub_command+="aprun -cc none -N 1 -d 23 -n $(( X * Y * Z )) ./espreso $EXAMPLEDIR/$EXAMPLE 0 ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"

   echo $qsub_command | tr ";" "\n" 
   echo $qsub_command | tr ";" "\n" | \
   qsub -q short -A $account -l select=$(( X * Y * Z )) -l walltime=00:10:00 -N $jobname
   #sbatch -N $(( X * Y * Z ))  -p test_large

  done

fi
