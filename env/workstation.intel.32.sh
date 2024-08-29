
. /opt/intel/oneapi/setvars.sh

. dependencies/install.gklib.sh icc
. dependencies/install.metis32.sh icc
. dependencies/install.parmetis32.sh mpiicc
. dependencies/install.precice.sh mpiicpc

export CXX=mpiicpc
