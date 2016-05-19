#!/bin/bash

if [ "$1" = "runpbs" ]; then



  export PARDISOLICMESSAGE=1
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs

#  export OMP_NUM_THREADS=1
#  export MKL_NUM_THREADS=1
#  export SOLVER_NUM_THREADS=24


export MKL_NUM_THREADS=23
export OMP_NUM_THREADS=23
export SOLVER_NUM_THREADS=23
export PAR_NUM_THREADS=23
export CILK_NWORKERS=23
export PARDISOLICMESSAGE=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs
export LC_CTYPE=""

  dom_size=$4
  clusters=$2
  clustt_size=$3
  corners=0

  d=${dom_size}
  c=${corners}

  x=${clustt_size}
  y=${clustt_size}
  z=${clustt_size}

  X=${clusters}
  Y=${clusters}
  Z=${clusters}

  log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c.log

  date | tee -a $log_file

  echo "Config: dom_size=  $d | cluster_size = $x:$y:$z | clusters = $X:$Y:$Z  "

  date | tee -a $log_file

  #mpirun -n $(( X * Y * Z )) ./espreso 0 ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file 

    cp -R ~/espreso/libs/        /scratch/work/user/lriha/esp/
    cp -R ~/espreso/examples/    /scratch/work/user/lriha/esp/
    cp    ~/espreso/espreso      /scratch/work/user/lriha/esp/
    cp    ~/espreso/salomon.sh   /scratch/work/user/lriha/esp/

    cd /scratch/work/user/lriha/esp/

 mpiexec_mpt -n $(( X * Y * Z )) perfboost -impi ./espreso examples/meshgenerator/cube_elasticity_fixed_bottom.txt 0 ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file

#  mpirun -f /home/lriha/espreso/hostfile perfboost -impi ./espreso examples/meshgenerator/cube_elasticity_fixed_bottom.txt 0 ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file

fi




#module load impi/5.0.3.048-iccifort-2015.3.187
#module load icc/2015.3.187
#module load imkl/11.2.3.187-iimpi-7.3.5

module load icc/2016.1.150-GCC-4.9.3 
module load imkl/11.3.1.150-iimpi-2016.00-GCC-4.9.3
module load tbb/4.4.2.152
module load DDT/5.0.1
module load itac/9.1.2.024

#INTEL MPI
#module load impi/5.1.2.150-iccifort-2016.1.150-GCC-4.9.3


#START - SGI MPI
module load perfboost
module load mpt/2.12

export MPI_DSM_DISTRIBUTE=0
export MPI_SHEPHERD=1
export PERFBOOST_VERBOSE=1
export MPI_VERBOSE=1
export MPI_BUFS_PER_PROC=512
#END -  SGI MPI



export MKL_NUM_THREADS=23
export OMP_NUM_THREADS=23
export SOLVER_NUM_THREADS=23
export PAR_NUM_THREADS=23
export CILK_NWORKERS=23
export PARDISOLICMESSAGE=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.
export LC_CTYPE=""

qsub_command_0="module load icc/2016.1.150-GCC-4.9.3;"
qsub_command_0+="module load imkl/11.3.1.150-iimpi-2016.00-GCC-4.9.3;"
qsub_command_0+="module load tbb/4.4.2.152;"
qsub_command_0+="module load DDT/5.0.1;"
qsub_command_0+="module load itac/9.1.2.024;"

#INTEL MPI
#qsub_command_0+="module load impi/5.1.2.150-iccifort-2016.1.150-GCC-4.9.3;"

#START - SGI MPI
qsub_command_0+="module load perfboost;"
qsub_command_0+="module load mpt/2.12;"

qsub_command_0+="export MPI_DSM_DISTRIBUTE=0;"
qsub_command_0+="export MPI_SHEPHERD=1;"
qsub_command_0+="export PERFBOOST_VERBOSE=1;"
qsub_command_0+="export MPI_VERBOSE=1;"
qsub_command_0+="export MPI_BUFS_PER_PROC=512;"
#END -  SGI MPI


qsub_command_0+="export MKL_NUM_THREADS=23;"
qsub_command_0+="export OMP_NUM_THREADS=23;"
qsub_command_0+="export SOLVER_NUM_THREADS=23;"
qsub_command_0+="export PAR_NUM_THREADS=23;"
qsub_command_0+="export CILK_NWORKERS=23;"
qsub_command_0+="export PARDISOLICMESSAGE=1;"
qsub_command_0+="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs:.;"
qsub_command_0+="export LC_CTYPE=;"


if [ "$#" -eq 0 ]; then
  echo "  Use one of the following commands:"
  echo "    ./salomon configure"
  echo "    ./salomon build"
  echo "    ./salomon run"
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
  ./waf configure --salomon
fi

if [ "$1" = "mesh" ]; then
  ./waf install --mesh
fi

if [ "$1" = "esmesh" ]; then
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs
  ./esmesh
fi

if [ "$1" = "build" ]; then
  ./waf install
fi

if [ "$1" = "build-pardiso_mkl" ]; then
  ./waf install --pardiso_mkl
fi

if [ "$1" = "mic" ]; then
  ./waf install --static --mic
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

  #               OM OK OK
  #               0   1   2   3   4   5   6   7   8   9   10
  dom_size=(      12  12  12  12  12  12  12  16  16  16  16)
  clustt_size_x=(  9   9   9   9   9   9   9  10  10  10  10)
# clustt_size_y=( 2   5   5   5   5   5   5   5   5   5 )
# clustt_size_z=( 1   5   5   5   5   5   5   5   5   5 )

  clusters_x=(    2   3   4   5   6   7   8   9   10  10  10)
  clusters_y=(    2   3   4   5   6   7   8   9    9  10  10)
  clusters_z=(    2   3   4   5   6   7   8   9    9   9  10)

  corners=(       0   0   0   0   0   0   0   0    0   0   0)

  for i in 0 1 2 3 4 5 6   # 8 9 10 8 9 10 # 0 1 2 3 4 5 6 7 8 9 10
  do
    d=${dom_size[${i}]}
    c=${corners[${i}]}

    x=${clustt_size_x[${i}]}
    y=${clustt_size_x[${i}]}
    z=${clustt_size_x[${i}]}

    X=${clusters_x[${i}]}
    Y=${clusters_y[${i}]}
    Z=${clusters_z[${i}]}

    jobname=espreso
    mpiranks=$(( X * Y * Z ))
    account=SERVICE
    actualTime=$( date +%y%m%d_%H:%M:%S )
    log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime.log
    out_dir=ESP-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c-$actualTime
    qsub_command=$qsub_command_0

    #date | tee -a $log_file
    qsub_command+="date | tee -a $log_file;"
    
    echo "Config: dom_size=  $d | cluster_size = $x:$y:$z | clusters = $X:$Y:$Z  "
    mkdir /scratch/work/user/lriha/esp/$out_dir
    cp -R ~/espreso/libs/     	 /scratch/work/user/lriha/esp/$out_dir
    cp -R ~/espreso/examples/ 	 /scratch/work/user/lriha/esp/$out_dir
    cp    ~/espreso/espreso      /scratch/work/user/lriha/esp/$out_dir
    cp    ~/espreso/salomon.sh   /scratch/work/user/lriha/esp/$out_dir
    cp    ~/espreso/salomon-qsub.sh /scratch/work/user/lriha/esp/$out_dir
    #cd /scratch/work/user/lriha/esp/


    #qsub_command+="mkdir /scratch/work/user/lriha/esp/$out_dir;"
    #qsub_command+="cp -R ~/espreso/libs/        /scratch/work/user/lriha/esp/$out_dir;"
    #qsub_command+="cp -R ~/espreso/examples/    /scratch/work/user/lriha/esp/$out_dir;"
    #qsub_command+="cp    ~/espreso/espreso      /scratch/work/user/lriha/esp/$out_dir;"
    #qsub_command+="cp    ~/espreso/salomon.sh   /scratch/work/user/lriha/esp/$out_dir;"
    #qsub_command+="cp    ~/espreso/salomon-qsub.sh /scratch/work/user/lriha/esp/$out_dir;"
    qsub_command+="cd /scratch/work/user/lriha/esp/$out_dir;"


    #mpirun -n $(( X * Y * Z )) hostname
    
   #INTEL MPI 
   #mpirun -n $(( X * Y * Z ))  ./espreso examples/meshgenerator/cube_elasticity_fixed_bottom.txt ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file
   #qsub_command+="mpirun -n $(( X * Y * Z ))  ./espreso examples/meshgenerator/cube_elasticity_fixed_bottom.txt ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file;"    

   #SGI MPT
   #mpiexec_mpt -n $mpiranks perfboost -impi ./espreso examples/meshgenerator/cube_elasticity_fixed_bottom.txt ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file
   qsub_command+="mpiexec_mpt -n $mpiranks perfboost -impi ./espreso examples/meshgenerator/cube_elasticity_fixed_bottom.txt ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file;"   
   echo $qsub_command
   echo $qsub_command | \
   qsub -A $account -q qmpp -l select=$mpiranks:ncpus=24:mpiprocs=1:ompthreads=1 -l walltime=00:20:00 -N $jobname

    #mpirun -n 10 ./espreso | tee -a $log_file
   
    #ddt -noqueue -start -n $(( X * Y * Z ))    ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} # | tee -a $log_file
   	
    #cp mesh.vtk mesh-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c.vtk
    #rm mesh.vtk
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
