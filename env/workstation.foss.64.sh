
buildname="build-workstation-foss-64"
export WAFLOCK=".lock-waf_linux_${buildname}"
rm -rf build
ln -s "${buildname}" build



. env/dependencies/install.suitesparse.sh foss gcc gfortran
. env/dependencies/install.gklib.sh foss gcc
. env/dependencies/install.metis64.sh foss gcc
. env/dependencies/install.parmetis64.sh foss mpicc
. env/dependencies/install.precice.sh foss mpic++

export CXX=mpic++
export ES_INT_WIDTH=64
