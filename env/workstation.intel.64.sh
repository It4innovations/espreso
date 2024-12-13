
. /opt/intel/oneapi/setvars.sh

. env/dependencies/install.gklib.sh icx ifx
. env/dependencies/install.metis64.sh icx
. env/dependencies/install.parmetis64.sh mpiicx
. env/dependencies/install.precice.sh mpiicpx

export CXX=mpiicpx
export ES_INT_WIDTH=64
