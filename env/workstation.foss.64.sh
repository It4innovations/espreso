
. env/dependencies/install.suitesparse.sh gcc gfortran
. env/dependencies/install.gklib.sh gcc
. env/dependencies/install.metis64.sh gcc
. env/dependencies/install.parmetis64.sh mpicc
. env/dependencies/install.precice.sh mpic++

export CXX=mpic++
export ES_INT_WIDTH=64
