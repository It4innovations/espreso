#!/bin/bash

WORKDIR=~/espreso/tests-jan2017/   #/scratch/temp/lriha/xeonphitemp-2  #ESCO2016-experiments/LargeTests/feti-static-3 #3-espreso-results-large-htfeti-FACEAVG-2clustersPerNode-lumped-11thredsPerMPI #FACEAVG
ESPRESODIR=~/espreso
EXAMPLEDIR=examples/meshgenerator
EXAMPLE=cube_elasticity_fixed_bottom.txt
THREADS_PER_MPI=11
MPI_PER_NODE=2
USE_MIC_NODES_ONLY=1
account=OPEN-7-46



#module load impi/5.0.3.048-iccifort-2015.3.187
#module load icc/2015.3.187
#module load imkl/11.2.3.187-iimpi-7.3.5

#module load icc/2016.1.150-GCC-4.9.3 
#module load imkl/11.3.1.150-iimpi-2016.00-GCC-4.9.3
#module load tbb/4.4.2.152
#module load DDT/5.0.1
#module load itac/9.1.2.024

#module load zlib/1.2.8
#module load CMake/3.5.2
#module load tbb/4.4.2.152
#module load intel/2017.00

module load impi/5.1.2.150-iccifort-2016.1.150-GCC-4.9.3-2.25
module load imkl/11.3.1.150-iimpi-2016.01-GCC-4.9.3-2.25 
#. /apps/all/imkl/2017.0.098-iimpi-2017.00-GCC-5.4.0-2.26/mkl/bin/mklvars.sh intel64
module load tbb/4.4.2.152 
module load zlib/1.2.8-intel-2016.01
module load CMake/3.3.1-intel-2016.01

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

if [ "$1" = "help" ]; then
  ./espreso -hh
fi

if [ "$1" = "configure" ]; then
  ./waf configure
  ./waf install
fi

if [ "$1" = "mesh" ]; then
  ./waf install -v
fi

if [ "$1" = "build" ]; then
  module list
  ./waf install 
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

  # *** Here you predefine the test you wantto execute *************
  #
  #               
  # test ID:      0   1   2   3   4   5   6   7   8   9
  dom_size=(      7   15  15  15  15  15  6  16  13  14 )
  clustt_size_x=( 6   11  11  11  11  11  5  10  10  10 )
  clustt_size_y=( 2   5   5   5   5   5   5   5   5   5 )
  clustt_size_z=( 1   5   5   5   5   5   5   5   5   5 )

  clusters_x=(    4   3   4   5   6   7   8  10   10   1 )
  clusters_y=(    4   3   4   5   6   7   8   9   10   1 )
  clusters_z=(    4   3   4   5   6   7   8   9   10   1 )

  #START - SGI MPI
  if [ "$2" = "sgi" ]; then
    module unload impi
    module load perfboost
    module load mpt/2.12

    export MPI_DSM_DISTRIBUTE=0
    export MPI_SHEPHERD=1
    export PERFBOOST_VERBOSE=1
    export MPI_VERBOSE=1
    export MPI_BUFS_PER_PROC=512
  fi
  #END -  SGI MPI

  qsub_command_0="#!/bin/bash;"
  qsub_command_0+="export MKL_NUM_THREADS=1;"
  qsub_command_0+="export OMP_NUM_THREADS=$THREADS_PER_MPI;"
  qsub_command_0+="export SOLVER_NUM_THREADS=$THREADS_PER_MPI;"
  qsub_command_0+="export PAR_NUM_THREADS=$THREADS_PER_MPI;"
  qsub_command_0+="export CILK_NWORKERS=1;"
  qsub_command_0+="export PARDISOLICMESSAGE=1;"
  qsub_command_0+="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.;"
  qsub_command_0+="export LC_CTYPE=;"
  qsub_command_0+="export MIC_ENV_PREFIX=MIC;"
  qsub_command_0+="export MIC_OMP_NUM_THREADS=60;"
  qsub_command_0+="export OFFLOAD_INIT=on_start;"
  qsub_command_0+="export MIC_USE_2MB_BUFFERS=10k;"
  qsub_command_0+="export MIC_OMP_NESTED=TRUE;"
  qsub_command_0+="export MIC_MKL_DYNAMIC=FALSE;"
  qsub_command_0+="export MIC_MKL_NUM_THREADS=3;"
  qsub_command_0+="export MIC_OMP_PROC_BIND=spread,close;"

  #if INTEL MPI
  if [ "$2" = "intel" ]; then
    qsub_command_0+="module load impi/5.1.2.150-iccifort-2016.1.150-GCC-4.9.3-2.25;"
  fi

  qsub_command_0+="module load imkl/11.3.1.150-iimpi-2016.01-GCC-4.9.3-2.25;"
  qsub_command_0+="module load tbb/4.4.2.152;"

  # if SGI MPI
  if [ "$2" = "sgi" ]; then
    qsub_command_0+="module unload impi;"
    qsub_command_0+="module unload iimpi;"

    qsub_command_0+="module load perfboost;"
    qsub_command_0+="module load mpt/2.12;"

    qsub_command_0+="export MPI_DSM_DISTRIBUTE=0;"
    qsub_command_0+="export MPI_SHEPHERD=1;"
    qsub_command_0+="export PERFBOOST_VERBOSE=1;"
    qsub_command_0+="export MPI_VERBOSE=1a;"
    qsub_command_0+="export MPI_BUFS_PER_PROC=512;"
  fi
  # end SGI MPI 
  
  qsub_command_0+="module list;"

  for i in 1 2 # 2 3 4 5 6 7 8 9 10
  do
    d=${dom_size[${i}]}

    x=${clustt_size_x[${i}]}
    y=${clustt_size_y[${i}]}
    z=${clustt_size_z[${i}]}

    X=${clusters_x[${i}]}
    Y=${clusters_y[${i}]}
    Z=${clusters_z[${i}]}


###  Overriding the table settings 

#    Scalability with clusters - changing the number of clusters/MPI processes 
    d=10 # subdomains size
    
    x=8 # cluster size in domains
    y=$x
    z=$x

    X=$i #$i
    Y=$X
    Z=$X
   
###  Cluster scalability - changing number of domains per cluster 
#    d=11 # subdomains size

#    x=$i # cluster size in domains
#    y=$x
#    z=$x

#    X=1
#    Y=$X
#    Z=$X

    # END overriding the 

    jobname=espreso
    mpiranks=$(( X * Y * Z ))
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
    cp       $ESPRESODIR/*.config       $WORKDIR/$out_dir

    qsub_command+="cd $WORKDIR/$out_dir;"
    qsub_command+="cat $¡PBS_NODEFILE | tee -a $node_file;"

    if [ "$3" = "mpi" ]; then

      export LD_LIBRARY_PATH=$ESPRESODIR/libs:$LD_LIBRARY_PATH  # /home/lriha/espreso/libs:$LD_LIBRARY_PATH
      export MKL_NUM_THREADS=1
      export OMP_NUM_THREADS=$THREADS_PER_MPI
      export SOLVER_NUM_THREADS=$THREADS_PER_MPI
      export PAR_NUM_THREADS=$THREADS_PER_MPI
      export CILK_NWORKERS=1
      export PARDISOLICMESSAGE=1
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.
      export MIC_ENV_PREFIX=MIC
      export MIC_OMP_NUM_THREADS=60
      export MIC_OMP_NESTED=TRUE
      export MIC_MKL_DYNAMIC=FALSE
      export MIC_OMP_PROC_BIND=spread,close
      export MIC_MKL_NUM_THREADS=3
      export MIC_USE_2MB_BUFFERS=100k
      export OFFLOAD_INIT=on_start
      export LC_CTYPE=

      cd $WORKDIR/$out_dir

      #cat $PBS_NODEFILE | tee -a $node_file
      if [ "$2" = "intel" ]; then
        mpirun      -n $(( X * Y * Z ))                ./espreso  $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file
      fi
      
      if [ "$2" = "sgi" ]; then
        mpiexec_mpt -n $(( X * Y * Z )) perfboost -impi  ./espreso $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file
      fi
    fi


    if [ "$3" = "pbs" ]; then
      cd $WORKDIR/$out_dir
      if [ "$2" = "sgi" ]; then
        log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime.warmup
        qsub_command+="mpiexec_mpt -n $(( X * Y * Z )) perfboost -impi ./espreso $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"
        log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime-run1.log
        qsub_command+="mpiexec_mpt -n $(( X * Y * Z )) perfboost -impi ./espreso $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"
        log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime-run2.log
        qsub_command+="mpiexec_mpt -n $(( X * Y * Z )) perfboost -impi ./espreso $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"
      fi
      if [ "$2" = "intel" ]; then
        log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime.warmup
        qsub_command+="mpirun      -n $(( X * Y * Z ))                 ./espreso $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"
        log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime-run1.log
        qsub_command+="mpirun      -n $(( X * Y * Z ))                 ./espreso $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"
        log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime-run2.log
        qsub_command+="mpirun      -n $(( X * Y * Z ))                 ./espreso $EXAMPLEDIR/$EXAMPLE ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"    
      fi

      echo $qsub_command | tr ";" "\n" | tr -d "¡" | tee $WORKDIR/$out_dir/job.qsub
      
      if [ "$USE_MIC_NODES_ONLY" = "1" ]; then
        #submit only to nodes with Xeon Phi 
        echo "qsub -q qmpp -A $account -l select=$(( ((X * Y * Z)/MPI_PER_NODE) + 1 )):ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI:accelerator=True -l walltime=00:20:00 -N $jobname" | tee $WORKDIR/$out_dir/submit_job_cmd.sh 
        echo $qsub_command | tr ";" "\n" | tr -d "¡" | \
        qsub -q qmpp -A $account -l select=$(( ((X * Y * Z)/MPI_PER_NODE) + 1 )):ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI:accelerator=True -l walltime=00:20:00 -N $jobname
      else 
        #submit to any nodes 
        echo "qsub -q qmpp -A $account -l select=$(( ((X * Y * Z)/MPI_PER_NODE) + 1 )):ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI -l walltime=00:20:00 -N $jobname" | tee $WORKDIR/$out_dir/submit_job_cmd.sh
        echo $qsub_command | tr ";" "\n" | tr -d "¡" | \
        qsub -q qmpp -A $account -l select=$(( ((X * Y * Z)/MPI_PER_NODE) + 1 )):ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI -l walltime=00:20:00 -N $jobname
      fi 
      
    fi

  done

fi



if [ "$1" = "run_ansys" ]; then

. ansys_run.conf

  qsub_command_0="#!/bin/bash;"
  qsub_command_0+="export MKL_NUM_THREADS=1;"
  qsub_command_0+="export OMP_NUM_THREADS=$THREADS_PER_MPI;"
  qsub_command_0+="export SOLVER_NUM_THREADS=$THREADS_PER_MPI;"
  qsub_command_0+="export PAR_NUM_THREADS=$THREADS_PER_MPI;"
  qsub_command_0+="export CILK_NWORKERS=1;"
  qsub_command_0+="export PARDISOLICMESSAGE=1;"
  qsub_command_0+="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.;"

  qsub_command_0+="export LC_CTYPE=;"

  qsub_command_0+="export MIC_ENV_PREFIX=MIC;"
  qsub_command_0+="export MIC_OMP_NUM_THREADS=60;"
  qsub_command_0+="export OFFLOAD_INIT=on_start;"
  qsub_command_0+="export MIC_USE_2MB_BUFFERS=10k;"
  qsub_command_0+="export MIC_OMP_NESTED=TRUE;"
  qsub_command_0+="export MIC_MKL_DYNAMIC=FALSE;"
  qsub_command_0+="export MIC_MKL_NUM_THREADS=3;"
  qsub_command_0+="export MIC_OMP_PROC_BIND=spread,close;"

  if [ "$2" = "sgi" ]; then
    #START - SGI MPI
    module unload impi
    module unload iimpi
    module load perfboost
    module load mpt/2.12

    export MPI_DSM_DISTRIBUTE=0
    export MPI_SHEPHERD=1
    export PERFBOOST_VERBOSE=1
    export MPI_VERBOSE=1
    export MPI_BUFS_PER_PROC=512
    #END -  SGI MPI
  fi

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

  # Select 1.) Number of nodes and number of rans per node 
#  for NODES in 32 #16 32 64 128 256 # compute nodes 
#  do
#   for MPI_PER_NODE in 22 #1 #2 6 12 24   # MPI processes per node 
#   do 
#    THREADS_PER_MPI=$(( (24)/(MPI_PER_NODE) ))    # number of threads per MPI process 
#    echo "THREADS_PER_MPI = " $THREADS_PER_MPI
#    RANKS=$(( (NODES)*(MPI_PER_NODE) ))      # MPI ranks

  # OR 
  # Select 2.) define number of ranks ... 
  #for RANKS in 1200 2400 4800 9600  #16 32 64 128 256 # compute nodes 
  for RANKS in "${MPIRANKS[@]}"
  do
   #for MPI_PER_NODE in 22
   #do
   # THREADS_PER_MPI=1
    NODES=$(( 1 + (RANKS)/(MPI_PER_NODE) ))
   
  # END 

    name=$EXAMPLE"_"$RANKS"_"$MPI_PER_NODE"_"$THREADS_PER_MPI
    echo "Directory name is: " $name

    actualTime=$( date +%y%m%d_%H:%M:%S )
    log_file=LOG-$name"_"$actualTime.log
    node_file=LOG-$name"_"$actualTime.node
    out_dir=ESP-$name"_"$actualTime
    qsub_command=$qsub_command_0
    qsub_command+="date | tee -a $log_file;"
    
    jobname=es_ans

    exadir=$EXAMPLE_DIR/$EXAMPLE"_"$RANKS
    echo $exadir

    mkdir  -p $WORKDIR
    mkdir  -p $WORKDIR/$out_dir
    cp -R    $ESPRESODIR/libs/               $WORKDIR/$out_dir
    cp       $ESPRESODIR/espreso             $WORKDIR/$out_dir
    cp       $ESPRESODIR/salomon.sh          $WORKDIR/$out_dir
    cp       $ESPRESODIR/*.config            $WORKDIR/$out_dir
    cp       $ESPRESODIR/ansys_run.conf      $WORKDIR/$out_dir

    qsub_command+="cd $WORKDIR/$out_dir;"
    qsub_command+="cat $¡PBS_NODEFILE | tee -a $node_file;"

    if [ "$3" = "mpi" ]; then

      # CPU Enviroment setup
      export MKL_NUM_THREADS=1
      export OMP_NUM_THREADS=$THREADS_PER_MPI
      export SOLVER_NUM_THREADS=$THREADS_PER_MPI
      export PAR_NUM_THREADS=$THREADS_PER_MPI
      export CILK_NWORKERS=1
      export PARDISOLICMESSAGE=1
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.

      # MIC environment setup 
      export MIC_ENV_PREFIX=MIC
      export MIC_OMP_NUM_THREADS=60
      export MIC_OMP_NESTED=TRUE
      export MIC_MKL_DYNAMIC=FALSE
      export MIC_KMP_AFFINITY=balanced
      export MIC_OMP_PROC_BIND=spread,close
      export MIC_MKL_NUM_THREADS=3
      export MIC_USE_2MB_BUFFERS=100k
      export OFFLOAD_INIT=on_start

      #OS X fix
      export LC_CTYPE=

      cd $WORKDIR/$out_dir
      cat $PBS_NODEFILE > $node_file

      if [ "$2" = "intel" ]; then
      	for f in "${FILES[@]}"
        do
          for DOM_PER_CLUST in "${SUBDOMAINS[@]}"
          do
            log_file=LOG-$name"_"$f"_"$DOM_PER_CLUST"_"$actualTime.log
            echo "mpirun -n $RANKS                       ./espreso -i esdata -p $exadir -c $f --SUBDOMAINS=$DOM_PER_CLUST               -vvv -mmm | tee -a $log_file"
	          mpirun -n $RANKS                             ./espreso -i esdata -p $exadir -c $f --SUBDOMAINS=$DOM_PER_CLUST               -vvv -mmm | tee -a $log_file
          done
	      done
      fi

      if [ "$2" = "sgi" ]; then
	      for f in "${FILES[@]}"
        do
          for DOM_PER_CLUST in "${SUBDOMAINS[@]}"
          do
            log_file=LOG-$name"_"$f"_"$DOM_PER_CLUST"_"$actualTime.log
            echo "mpiexec_mpt -n $RANKS perfboost -impi  ./espreso -i esdata -p $exadir -c $f --SUBDOMAINS=$DOM_PER_CLUST               -vvv -mmm | tee -a $log_file"
                  mpiexec_mpt -n $RANKS perfboost -impi  ./espreso -i esdata -p $exadir -c $f --SUBDOMAINS=$DOM_PER_CLUST               -vvv -mmm | tee -a $log_file
          done
        done
      fi
    
    fi # end of MPI run 


    if [ "$3" = "pbs" ]; then
      cd $WORKDIR/$out_dir

      if [ "$2" = "sgi" ]; then
        f=${FILES[0]}
        DOM_PER_CLUST=${SUBDOMAINS[0]}
        log_file=LOG-$name"_"$f"_"$DOM_PER_CLUST"_"$actualTime.warmup
        qsub_command+="mpiexec_mpt -n $RANKS perfboost -impi ./espreso -i esdata -p ${exadir} -c ${f} --SUBDOMAINS=$DOM_PER_CLUST  -vvv -mmm | tee -a $log_file;"

	      for f in "${FILES[@]}"
        do
          for DOM_PER_CLUST in "${SUBDOMAINS[@]}"
          do 
            log_file=LOG-$name"_"$f"_"$DOM_PER_CLUST"_"$actualTime-run1.log
            qsub_command+="mpiexec_mpt -n $RANKS perfboost -impi ./espreso -i esdata -p ${exadir} -c ${f} --SUBDOMAINS=$DOM_PER_CLUST  -vvv -mmm | tee -a $log_file;"
            log_file=LOG-$name"_"$f"_"$DOM_PER_CLUST"_"$actualTime-run2.log
            qsub_command+="mpiexec_mpt -n $RANKS perfboost -impi ./espreso -i esdata -p ${exadir} -c ${f} --SUBDOMAINS=$DOM_PER_CLUST  -vvv -mmm | tee -a $log_file;"
	        done
        done
      fi

      if [ "$2" = "intel" ]; then
        f=${FILES[0]}
        DOM_PER_CLUST=${SUBDOMAINS[0]}
        log_file=LOG-$name"_"$f"_"$DOM_PER_CLUST"_"$actualTime.warmup
        qsub_command+="mpirun      -n $RANKS                 ./espreso -i esdata -p ${exadir} -c ${f} --SUBDOMAINS=$DOM_PER_CLUST  -vvv -mmm | tee -a $log_file;"

	      for f in "${FILES[@]}"
	      do
          for DOM_PER_CLUST in "${SUBDOMAINS[@]}"
          do
            log_file=LOG-$name"_"$f"_"$DOM_PER_CLUST"_"$actualTime-run1.log
            qsub_command+="mpirun      -n $RANKS                 ./espreso -i esdata -p ${exadir} -c ${f} --SUBDOMAINS=$DOM_PER_CLUST  -vvv -mmm | tee -a $log_file;"
            log_file=LOG-$name"_"$f"_"$DOM_PER_CLUST"_"$actualTime-run2.log
            qsub_command+="mpirun      -n $RANKS                 ./espreso -i esdata -p ${exadir} -c ${f} --SUBDOMAINS=$DOM_PER_CLUST  -vvv -mmm | tee -a $log_file;"
          done
	      done
      fi

      echo $qsub_command | tr ";" "\n" | tr -d "¡" | tee $WORKDIR/$out_dir/job.qsub

      if [ "$USE_MIC_NODES_ONLY" = "1" ]; then
        #submit only to nodes with Xeon Phi 
        echo "qsub -q $QUEUE -A $account -l select=$NODES:ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI:accelerator=True -l walltime=00:60:00 -N $jobname" | tee $WORKDIR/$out_dir/submit_job_cmd.sh
        echo $qsub_command | tr ";" "\n" | tr -d "¡" | \
            qsub -q $QUEUE -A $account -l select=$NODES:ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI:accelerator=True  -l walltime=00:60:00 -N $jobname

      else 
        #submit to any type of nodes 
        echo "qsub -q $QUEUE -A $account -l select=$NODES:ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI -l walltime=00:60:00 -N $jobname" | tee $WORKDIR/$out_dir/submit_job_cmd.sh
        echo $qsub_command | tr ";" "\n" | tr -d "¡" | \
            qsub -q $QUEUE -A $account -l select=$NODES:ncpus=24:mpiprocs=$MPI_PER_NODE:ompthreads=$THREADS_PER_MPI -l walltime=00:60:00 -N $jobname

      fi 

    fi

   #done 
  done

fi








if [ "$1" = "debug" ]; then
  module load valgrind/3.9.0-impi
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs
  export OMP_NUM_THREADS=1
  export MKL_PARDISO_OOC_MAX_CORE_SIZE=3500
  export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000

  for i in 0 #1 2 3 do
  do

    log_file=LOG-1:1:1-2:2:2-5:5:5.log

    date | tee $log_file

    echo "Config: dom_size = 5 | cluster_size = 2:2:2 | clusters = 1:1:1"

    date | tee -a $log_file

    valgrind ./espreso ${el_type[${i}]} 1 1 1 2 2 2 5 5 5 | tee -a $log_file
  done
fi
