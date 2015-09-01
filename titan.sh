#!/bin/bash

module swap PrgEnv-pgi/5.2.40 PrgEnv-intel/5.2.40
#module unload cray-libsci/13.0.4
module swap intel/14.0.2.144  intel/15.0.2.164 #  intel/13.1.3.192 # intel/15.0.2.164
module unload cray-libsci
module load gcc
module load cray-tpsl/1.5.0
#module load metis/5.1.0
module list

#cc -craype-verbose
#CC -craype-verbose
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../libs


#cd build
#CC -Wall -openmp mesh/src/settings.cpp.1.o mesh/src/loader.cpp.1.o mesh/src/structures/coordinates.cpp.1.o mesh/src/structures/mesh.cpp.1.o mesh/src/structures/boundaries.cpp.1.o mesh/src/matrices/matrix.cpp.1.o mesh/src/matrices/denseMatrix.cpp.1.o mesh/src/elements/element.cpp.1.o mesh/src/elements/3D/hexahedron8.cpp.1.o mesh/src/elements/3D/hexahedron20.cpp.1.o mesh/src/elements/3D/tetrahedron4.cpp.1.o mesh/src/elements/3D/tetrahedron10.cpp.1.o mesh/src/elements/3D/prisma6.cpp.1.o mesh/src/elements/3D/prisma15.cpp.1.o mesh/src/elements/3D/pyramid5.cpp.1.o mesh/src/elements/3D/pyramid13.cpp.1.o mesh/src/elements/2D/square.cpp.1.o mesh/src/elements/2D/triangle.cpp.1.o mesh/src/elements/1D/line.cpp.1.o mesh/src/elements/1D/point3d.cpp.1.o mesh/src/elements/1D/point2d.cpp.1.o mesh/src/main.cpp.1.o -o /autofs/nccs-svm1_home1/lriha/espreso-lriha/build/mesh/esmesh -Wl,-Bdynamic -L../libs -Lbem -lmetis32 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lesbem

#./mesh/esmesh
#exit

#export CPATH=/apps/icc/2015.3.187-GNU-5.1.0-2.25/tbb/include/:$CPATH

#module load impi/5.0.3.048-iccifort-2015.3.187
#module load MPT/2.12

#module load OpenMPI/1.8.6-GNU-5.1.0-2.25
#module unload GCC/5.1.0-binutils-2.25 binutils/2.25-GCC-5.1.0-binutils-2.25 GNU/5.1.0-2.25 impi/5.0.3.048-iccifort-2015.3.187

#module load icc/2015.3.187
#module load imkl/11.2.3.187-iimpi-7.3.5
#module load DDT/5.0.1

#export OMPI_CXX=icpc

#mpic++ -V

#module load CMake/3.0.0-intel-2015b 
#module load tbb/4.3.5.187
export LC_CTYPE=""

if [ "$#" -ne 1 ]; then
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
  ./waf configure --titan --static
fi

if [ "$1" = "build" ]; then
  ./waf install -v
#cd tools/metis-5.1.0/ ; make config compiler=cc #; make
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
  export PARDISOLICMESSAGE=1 
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs

  export OMP_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export SOLVER_NUM_THREADS=24

  export MKL_PARDISO_OOC_MAX_CORE_SIZE=3500
  export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs/

  log_file=output_8_nodes.log

  dom_size=10      #20
  clusters=3       #2
  domains_p_c=5    #7



  #               OM OK OK
  #               0   1   2   3   4   5   6   7   8   9
  dom_size=(      16  16  16  16  16  16  16  16  13  14 )
  clustt_size_x=( 10  10  10  10  10  10  10  10  10  10 )
# clustt_size_y=( 2   5   5   5   5   5   5   5   5   5 )
# clustt_size_z=( 1   5   5   5   5   5   5   5   5   5 )

  clusters_x=(    2   3   4   5   6   7   8   9   1   1 )
  clusters_y=(    2   3   4   5   6   7   8   9   1   1 )
  clusters_z=(    2   3   4   5   6   7   8   9   1   1 )

  corners=(       0   0   0   0   0   0   0   0   0   0 )

  for i in 7 # 0 1 2 3 4 5 6
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

    date | tee $log_file

    echo "Config: dom_size=  $d | cluster_size = $x:$y:$z | clusters = $X:$Y:$Z  "

    date | tee -a $log_file

    #                                          ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}   | tee -a $log_file
    #mpirun -bind-to-none -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}   | tee -a $log_file
    #mpirun -bind-to none -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}               # | tee -a $log_file
   
    aprun -n 1 ./esmesh
    #mpirun -n $(( X * Y * Z ))  ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d}                   | tee -a $log_file

   
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
