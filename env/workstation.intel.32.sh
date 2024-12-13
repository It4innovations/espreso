
. /opt/intel/oneapi/setvars.sh

. env/dependencies/install.gklib.sh icx ifx
. env/dependencies/install.metis32.sh icx
. env/dependencies/install.parmetis32.sh mpiicx
. env/dependencies/install.precice.sh mpiicpx

export CXX=mpiicpx
export ES_INT_WIDTH=32
