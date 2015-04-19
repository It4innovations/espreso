#!/bin/bash
##module load intel/13.5.192
#module load intel/14.0.1
#module load intel/15.2.164
#module load impi/4.1.1.036
#module load impi/4.1.0.030
#module load openmpi/1.8.1-icc

#. /apps/intel/BETA2015/haswell.sh

# Parameters :
#   - TM_BLOCK_START (use -DTM_BLOCK_START) - adds a MPI_Barrier before time measurement - barrier time is not measured
#   - TM_BLOCK_END   (use -DTM_BLOCK_END)   - adds a MPI_Barrier before end of time measurement - barrier time is part of the measured time
#   - USE_MPI_3      (use -DUSE_MPI_3)      - use non-blocking MPI_Iallreduce - in Pipeline CG algorithm - requires MPI 3 standard implementation

rm -fr bin
mkdir bin

cp Makefile.Anselm.devel bin/Makefile

find src -type f -name "*.cpp" -exec cp {} bin  \;
find src -type f -name "*.h" -exec cp {} bin  \;

cd bin
make cleanall
make lib -j 16
cd ..

#make libcuda -j 16
#cp libsolver_cuda.so libsolver.so

#make clean

#make clean
#make debug -j 16

#make cleanall
#make libmic -j 16
#make clean

#cp libsolver_mic.so ../gen/

