#!/bin/bash

WORKDIR=~/espreso-results-pbs-static-pbs
ESPRESODIR=~/espreso_git/espreso
EXAMPLEDIR=examples/meshgenerator
EXAMPLE=cube_elasticity_fixed_bottom.txt
THREADS_PER_MPI=24
MPI_PER_NODE=1

#module load impi/5.0.3.048-iccifort-2015.3.187
#module load icc/2015.3.187
#module load imkl/11.2.3.187-iimpi-7.3.5

#module load icc/2016.1.150-GCC-4.9.3 
#module load imkl/11.3.1.150-iimpi-2016.00-GCC-4.9.3
#module load tbb/4.4.2.152
#module load DDT/5.0.1
#module load itac/9.1.2.024

module load impi/5.1.2.150-iccifort-2016.1.150-GCC-4.9.3-2.25
module load imkl/11.3.1.150-iimpi-2016.01-GCC-4.9.3-2.25 
module load tbb/4.4.2.152 

if [ "$#" -eq 0 ]; then
  echo "  Use one of the following commands:"
  echo "    ./salomon configure"
  echo "    ./salomon build"
  echo "    ./salomon run sgi|intel mpi|pbs"
  echo "    ./salomon clean"
  echo "    ./salomon distclean"
  echo ""
  echo "  'configure' sets all parameters for the compiler. It is mandatory to run this command at the first time."
  echo "  'build' makes the runable application."
  echo "  'run' starts the computation. The result is put to Paraview input file 'mesh.vtk'."
  echo "  'clean' removes all files created by build process."
  echo "  'distclean' removes all files."
fi

if [ "$1" = "configure" ]; then
  ./waf configure
fi

if [ "$1" = "mesh" ]; then
  ./waf install -v
fi

if [ "$1" = "build" ]; then
  module list
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

  if [ "$2" = "sgi" ]; then
    #START - SGI MPI
    module unload impi
    module load perfboost
    module load mpt/2.12

    export MPI_DSM_DISTRIBUTE=0
    export MPI_SHEPHERD=1
    export PERFBOOST_VERBOSE=1
    export MPI_VERBOSE=1
    export MPI_BUFS_PER_PROC=512
    #END -  SGI MPI
  fi

  #               OM OK OK
  #               0   1   2   3   4   5   6   7   8   9
  dom_size=(      7   15  15  15  15  15  6  16  13  14 )
  clustt_size_x=( 6   11  11  11  11  11  5  10  10  10 )
  clustt_size_y=( 2   5   5   5   5   5   5   5   5   5 )
  clustt_size_z=( 1   5   5   5   5   5   5   5   5   5 )

  clusters_x=(    4   3   4   5   6   7   8  10   10   1 )
  clusters_y=(    4   3   4   5   6   7   8   9   10   1 )
  clusters_z=(    4   3   4   5   6   7   8   9   10   1 )

  corners=(       0   0   0   0   0   0   0   0   0   0 )

  qsub_command_0="#!/bin/bash;"
  qsub_command_0+="export MKL_NUM_THREADS=1;"
  qsub_command_0+="export OMP_NUM_THREADS=1;"
  qsub_command_0+="export SOLVER_NUM_THREADS=$THREADS_PER_MPI;"
  qsub_command_0+="export PAR_NUM_THREADS=$THREADS_PER_MPI;"
  qsub_command_0+="export CILK_NWORKERS=$THREADS_PER_MPI;"
  qsub_command_0+="export PARDISOLICMESSAGE=1;"
  qsub_command_0+="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.;"
  qsub_command_0+="export LC_CTYPE=;"

  if [ "$2" = "intel" ]; then
    qsub_command_0+="module load impi/5.1.2.150-iccifort-2016.1.150-GCC-4.9.3-2.25;"
  fi

  qsub_command_0+="module load imkl/11.3.1.150-iimpi-2016.01-GCC-4.9.3-2.25;"
  qsub_command_0+="module load tbb/4.4.2.152;"

  if [ "$2" = "sgi" ]; then
    qsub_command_0+="module unload impi;"
    qsub_command_0+="module unload iimpi;"

    qsub_command_0+="module load perfboost;"
    qsub_command_0+="module load mpt/2.12;"

    qsub_command_0+="export MPI_DSM_DISTRIBUTE=0;"
    qsub_command_0+="export MPI_SHEPHERD=1;"
    qsub_command_0+="export PERFBOOST_VERBOSE=1;"
    qsub_command_0+="export MPI_VERBOSE=1;"
    qsub_command_0+="export MPI_BUFS_PER_PROC=512;"
  fi
  
  qsub_command_0+="module list;"
	
  for i in 1 # 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
  do
    d=${dom_size[${i}]}
    c=${corners[${i}]}

    x=${clustt_size_x[${i}]}
    y=${clustt_size_y[${i}]}
    z=${clustt_size_z[${i}]}

    X=${clusters_x[${i}]}
    Y=${clusters_y[${i}]}
    Z=${clusters_z[${i}]}

    d=7 # subdomains size
    c=0 # number of corners - nefunguje 
    x=6 # cluster size in domains

    y=$x
    z=$x

    X=$i
    Y=$X
    Z=$X

    jobname=espreso
    mpiranks=$(( X * Y * Z ))
    account=SERVICE
    actualTime=$( date +%y%m%d_%H:%M:%S )
    log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime.log
    node_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime.node
    out_dir=ESP-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime
    qsub_command=$qsub_command_0
    qsub_command+="date | tee -a $log_file;"

    echo "Config: dom_size=  $d | cluster_size = $x:$y:$z | clusters = $X:$Y:$Z  "
    mkdir    $WORKDIR/$out_dir
    mkdir -p $WORKDIR/$out_dir/$EXAMPLEDIR
    cp -R    $ESPRESODIR/libs/          $WORKDIR/$out_dir
    cp -R    $ESPRESODIR/$EXAMPLEDIR/*  $WORKDIR/$out_dir/$EXAMPLEDIR
    cp       $ESPRESODIR/espreso        $WORKDIR/$out_dir
    cp       $ESPRESODIR/salomon.sh     $WORKDIR/$out_dir

    qsub_command+="cd $WORKDIR/$out_dir;"
    qsub_command+="cat $¡PBS_NODEFILE | tee -a $node_file;"

    if [ "$3" = "mpi" ]; then

      export MKL_NUM_THREADS=1
      export OMP_NUM_THREADS=1
      export SOLVER_NUM_THREADS=$THREADS_PER_MPI
      export PAR_NUM_THREADS=$THREADS_PER_MPI
      export CILK_NWORKERS=$THREADS_PER_MPI
      export PARDISOLICMESSAGE=1
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.
      export LC_CTYPE=

      cd $WORKDIR/$out_dir
      cat $PBS_NODEFILE | tee -a $node_file

      if [ "$2" = "intel" ]; then
        mpirun      -n $(( X * Y * Z ))                  ./espreso examples/meshgenerator/cube_elasticity_fixed_bottom.txt ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file
      fi
      
      if [ "$2" = "sgi" ]; then
        mpiexec_mpt -n $(( X * Y * Z )) perfboost -impi  ./espreso examples/meshgenerator/cube_elasticity_fixed_bottom.txt ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file
      fi
    fi


    if [ "$3" = "pbs" ]; then
      cd $WORKDIR/$out_dir
      if [ "$2" = "sgi" ]; then
        qsub_command+="mpiexec_mpt -n $(( X * Y * Z )) perfboost -impi ./espreso $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"
      fi
      if [ "$2" = "intel" ]; then
        qsub_command+="mpirun      -n $(( X * Y * Z ))                 ./espreso $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"    
      fi

      echo $qsub_command | tr ";" "\n" | tr -d "¡" | tee $WORKDIR/$out_dir/job.qsub
      echo "qsub -q qmpp -A $account -l select=$(( ((X * Y * Z)/MPI_PER_NODE) + 1 )):ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI -l walltime=00:20:00 -N $jobname" | tee -a $WORKDIR/$out_dir/job.qsub

      echo $qsub_command | tr ";" "\n" | tr -d "¡" | \
      qsub -q qmpp -A $account -l select=$(( ((X * Y * Z)/MPI_PER_NODE) + 1 )):ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI -l walltime=00:20:00 -N $jobname
    
    fi

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
