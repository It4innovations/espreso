#!/bin/bash

ml LUMI/24.03
ml buildtools/24.03
ml Boost/1.83.0-cpeCray-24.03
ml Eigen/3.4.0
ml buildtools-python/24.03-cray-python3.11

export CXX=CC
export ES_INT_WIDTH=32
export BLAS_LIBRARIES=sci_cray
export LAPACK_LIBRARIES=sci_cray

. dependencies/install.gklib.sh cc
. dependencies/install.metis32.sh cc
. dependencies/install.parmetis32.sh mpicc
. dependencies/install.suitesparse.sh cc "-DLAPACK_LIBRARIES=/opt/cray/pe/libsci/24.03.0/CRAYCLANG/17.0/x86_64/lib/libsci_cray.so -DBLAS_LIBRARIES=/opt/cray/pe/libsci/24.03.0/CRAYCLANG/17.0/x86_64/lib/libsci_cray.so"
. dependencies/install.precice.sh mpic++
