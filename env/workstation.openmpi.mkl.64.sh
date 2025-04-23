
buildname="build-workstation-openmpi-mkl-64"
export WAFLOCK=".lock-waf_linux_${buildname}"
rm -rf build
ln -s "${buildname}" build



. env/dependencies/install.suitesparse.sh mkl gcc gfortran
. env/dependencies/install.gklib.sh mkl gcc
. env/dependencies/install.metis64.sh mkl gcc
. env/dependencies/install.parmetis64.sh mkl mpicc
. env/dependencies/install.precice.sh mkl mpic++

export CXX=mpic++
export ES_INT_WIDTH=64
export CPATH=/usr/include/mkl:$CPATH
