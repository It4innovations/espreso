
buildname="build-workstation-intel-64"
export WAFLOCK=".lock-waf_linux_${buildname}"
if [ -d build ]; then
  curr_buildname="$(readlink -f build)"
fi
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi



. /opt/intel/oneapi/setvars.sh

. env/dependencies/install.gklib.sh intel icx ifx
. env/dependencies/install.metis64.sh intel icx
. env/dependencies/install.parmetis64.sh intel mpiicx
. env/dependencies/install.precice.sh intel mpiicpx

export CXX=mpiicpx
export ES_INT_WIDTH=64
