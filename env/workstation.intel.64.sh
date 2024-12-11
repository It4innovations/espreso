
. /opt/intel/oneapi/setvars.sh

. dependencies/install.gklib.sh icx
. dependencies/install.metis64.sh icx
. dependencies/install.parmetis64.sh mpiicx
. dependencies/install.precice.sh mpiicpx

export CXX=mpiicpx
export ES_INT_WIDTH=64
