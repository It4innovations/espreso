
buildname="build-workstation-foss-32"
export WAFLOCK=".lock-waf_linux_${buildname}"
rm -rf build
ln -s "${buildname}" build



. env/dependencies/install.suitesparse.sh foss gcc gfortran
. env/dependencies/install.gklib.sh foss gcc
. env/dependencies/install.metis32.sh foss gcc
. env/dependencies/install.parmetis32.sh foss mpicc
. env/dependencies/install.precice.sh foss mpic++

export CXX=mpic++
export ES_INT_WIDTH=32
