
. /opt/intel/oneapi/setvars.sh

. dependencies/install.metis32.sh icc
. dependencies/install.parmetis32.sh mpiicc

export CXX=mpiicpc
