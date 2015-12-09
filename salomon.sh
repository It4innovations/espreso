#!/bin/bash

#INTEL=/apps/all/icc/2015.3.187
#. $INTEL/tbb/bin/tbbvars.sh intel64

if [ "$1" = "runpbs" ]; then



  export PARDISOLICMESSAGE=1
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs

  export OMP_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export SOLVER_NUM_THREADS=24

  export MKL_PARDISO_OOC_MAX_CORE_SIZE=3500
  export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs/

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

  mpirun -f /home/lriha/espreso/hostfile perfboost -impi ./espreso 0 ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file

fi




#module load impi/5.0.3.048-iccifort-2015.3.187
#module load icc/2015.3.187
#module load imkl/11.2.3.187-iimpi-7.3.5

module load impi/5.1.2.150-iccifort-2016.1.150-GCC-4.9.3
module load icc/2016.1.150-GCC-4.9.3 
module load imkl/11.3.1.150-iimpi-2016.00-GCC-4.9.3
module load tbb/4.4.2.152
module load DDT/5.0.1
module load itac/9.1.2.024

export MKL_NUM_THREADS=12
export OMP_NUM_THREADS=12
export SOLVER_NUM_THREADS=12
export PAR_NUM_THREADS=12
export CILK_NWORKERS=12
export PARDISOLICMESSAGE=1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs
export LC_CTYPE=""

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
  #               0   1   2   3   4   5   6   7   8   9
  dom_size=(      11  11  13  16  16  16  16  13  13  14 )
  clustt_size_x=( 12  12  10  10  10  10  10  10  10  10 )
# clustt_size_y=( 2   5   5   5   5   5   5   5   5   5 )
# clustt_size_z=( 1   5   5   5   5   5   5   5   5   5 )

  clusters_x=(    2   3   4   5   6   7   8   9   10   1 )
  clusters_y=(    2   3   4   5   6   7   8   9   10   1 )
  clusters_z=(    2   3   4   5   6   7   8   9   10   1 )

  corners=(       0   0   0   0   0   0   0   0   0   0 )

  for i in 8 # 0 1 2 3 4 5 6
  do
    d=${dom_size[${i}]}
    c=${corners[${i}]}

    x=${clustt_size_x[${i}]}
    y=${clustt_size_x[${i}]}
    z=${clustt_size_x[${i}]}

    X=${clusters_x[${i}]}
    Y=${clusters_y[${i}]}
    Z=${clusters_z[${i}]}

    log_file=LOG-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c.log

    date | tee -a $log_file

    echo "Config: dom_size=  $d | cluster_size = $x:$y:$z | clusters = $X:$Y:$Z  "

    date | tee -a $log_file

    #                                          ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}   | tee -a $log_file
    #mpirun -bind-to-none -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}   | tee -a $log_file
    #mpirun -bind-to none -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}               # | tee -a $log_file
    
    cp -R ~/espreso/libs/     	 /scratch/work/user/lriha/esp/
    cp -R ~/espreso/examples/ 	 /scratch/work/user/lriha/esp/
    cp    ~/espreso/espreso      /scratch/work/user/lriha/esp/
    cp    ~/espreso/salomon.sh   /scratch/work/user/lriha/esp/
    
    cd /scratch/work/user/lriha/esp/
    mpirun -n $(( X * Y * Z )) hostname
    mpirun -n $(( X * Y * Z ))  ./espreso examples/meshgenerator/cube_elasticity_fixed_bottom.txt ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file
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
