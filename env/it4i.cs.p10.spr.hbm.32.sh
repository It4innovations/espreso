
buildname="build-it4i-cs-p10-spr-hbm-32"
export WAFLOCK=".lock-waf_linux_${buildname}"
curr_buildname="$(readlink build)"
if [ "${curr_buildname}" != "${buildname}" ]
then
    rm -rf build
    ln -s "${buildname}" build
fi

ml Python/3.10.4-GCCcore-11.3.0
ml OpenMPI/4.1.4-GCC-11.3.0
ml OpenBLAS/0.3.20-GCC-11.3.0
ml CMake/3.23.1-GCCcore-11.3.0

. env/dependencies/install.suitesparse.sh p10 gcc gfortran
. env/dependencies/install.gklib.sh p10 gcc
. env/dependencies/install.metis32.sh p10 gcc
. env/dependencies/install.parmetis32.sh p10 mpicc
. env/dependencies/install.nose2.sh p10

export CXX=mpic++
export ES_INT_WIDTH=32
export BLAS_LIBRARIES=openblas
export LAPACK_LIBRARIES=openblas
