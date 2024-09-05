
ml intel/2023b CMake/3.27.6-GCCcore-13.2.0 Boost/1.83.0-GCC-13.2.0 Eigen/3.4.0-GCCcore-13.2.0 libxml2/2.11.5-GCCcore-13.2.0 Python/3.11.5-GCCcore-13.2.0 Boost.Python-NumPy/1.83.0-gfbf-2023b

. dependencies/install.suitesparse.sh icx
. dependencies/install.gklib.sh icx
. dependencies/install.metis32.sh icx
. dependencies/install.parmetis32.sh mpiicx
. dependencies/install.precice.sh mpiicpx

export CXX=mpiicpx
