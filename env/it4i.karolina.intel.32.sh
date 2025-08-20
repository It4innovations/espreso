
buildname="build-it4i-karolina-intel-32"
export WAFLOCK=".lock-waf_linux_${buildname}"
curr_buildname="$(readlink build)"
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi



ml intel/2024a CMake/3.27.6-GCCcore-13.2.0 Boost/1.83.0-GCC-13.2.0 Eigen/3.4.0-GCCcore-13.2.0 libxml2/2.11.5-GCCcore-13.2.0 Python/3.11.5-GCCcore-13.2.0 Boost.Python-NumPy/1.83.0-gfbf-2023b

. env/dependencies/install.suitesparse.sh intel icx ifx
. env/dependencies/install.gklib.sh intel icx
. env/dependencies/install.metis32.sh intel icx
. env/dependencies/install.parmetis32.sh intel mpiicx
. env/dependencies/install.precice.sh intel mpiicpx

export CXX=mpiicpx
export ES_INT_WIDTH=32
