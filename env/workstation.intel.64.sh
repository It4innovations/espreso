
. /opt/intel/oneapi/setvars.sh

. dependencies/install.gklib.sh icc
. dependencies/install.metis64.sh icc
. dependencies/install.parmetis64.sh mpiicc
. dependencies/install.precice.sh mpiicpc

export CXX=mpiicpc