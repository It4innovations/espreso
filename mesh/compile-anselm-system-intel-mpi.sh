#!/bin/bash
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

cp Makefile.Anselm bin/Makefile

cp -R src/* bin

cp ../solver/bin/libsolver.so libs

cd bin
make cleanall
make opt -j 16
#make debug -j 16
make clean
cd ..

cp bin/esmesh ../libs
#cp pmcube_db ../


#make cleanall
#make mic -j 16
#make clean

#cp pmcube_mic ../

#rm *.cpp
#rm *.h
