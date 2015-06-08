#!/bin/bash

module load bullxmpi/bullxmpi_1.2.4.1
#module load impi/4.1.1.036
module load intel/15.2.164
module load cmake/2.8.11
export LC_CTYPE=""

if [ "$#" -ne 1 ]; then
  echo "  Use one of the following commands:"
  echo "    ./anselm configure"
  echo "    ./anselm build"
  echo "    ./anselm run"
  echo "    ./anselm clean"
  echo "    ./anselm distclean"
  echo ""
  echo "  'configure' sets all parameters for the compiler. It is mandatory to run this command at the first time."
  echo "  'build' makes the runable application."
  echo "  'run' starts the computation. The result is put to Paraview input file 'mesh.vtk'."
  echo "  'clean' removes all files created by build process."
  echo "  'distclean' removes all files."
fi

if [ "$1" = "configure" ]; then
  ./waf configure --anselm
fi

if [ "$1" = "build" ]; then
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

if [ "$1" = "run" ]; then
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./libs
  export OMP_NUM_THREADS=1
  export MKL_PARDISO_OOC_MAX_CORE_SIZE=3500
  export MKL_PARDISO_OOC_MAX_SWAP_SIZE=2000

  log_file=output_8_nodes.log

  dom_size=10      #20
  clusters=3       #2
  domains_p_c=5    #7

  #         HEXA8 TETRA4 TETRA10
  el_type=(   0     1      2    )

  #               OM OK OK
  #               0   1   2   3   4   5   6   7   8   9
  dom_size=(      15  16  17  8   9  10  11  12  13  14 )
  clustt_size_x=( 10  10  10  10  10 10  10  10  10  10 )

  clustt_size_y=( 5   5   5   5   5   5   5   5   5   5 )
  clustt_size_z=( 5   5   5   5   5   5   5   5   5   5 )

  clusters_x=(    1   1   1   1   1   1   1   1   1   1 )
  clusters_y=(    1   1   1   1   1   1   1   1   1   1 )
  clusters_z=(    1   1   1   1   1   1   1   1   1   1 )

  corners=(       0   0   0   0   0   0   0   0   0   0 )

  for i in 1 2 # 0 1 2 3 4 5 6
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

    ./espreso ${el_type[0]} ${X} ${Y} ${Z} ${x} ${y} ${z} ${d} ${d} ${d} | tee -a $log_file

    cp mesh.vtk mesh-$X:$Y:$Z-$x:$y:$z-$d:$d:$d-$c:$c:$c.vtk
    rm mesh.vtk
  done

fi
