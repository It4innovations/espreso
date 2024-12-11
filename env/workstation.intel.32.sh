
. /opt/intel/oneapi/setvars.sh

. dependencies/install.gklib.sh icx
. dependencies/install.metis32.sh icx
. dependencies/install.parmetis32.sh mpiicx
. dependencies/install.precice.sh mpiicpx

export CXX=mpiicpx
export ES_INT_WIDTH=32
