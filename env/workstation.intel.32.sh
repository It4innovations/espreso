
. /opt/intel/oneapi/setvars.sh

. env/dependencies/install.gklib.sh intel icx ifx
. env/dependencies/install.metis32.sh intel icx
. env/dependencies/install.parmetis32.sh intel mpiicx
. env/dependencies/install.precice.sh intel mpiicpx

export CXX=mpiicpx
export ES_INT_WIDTH=32
