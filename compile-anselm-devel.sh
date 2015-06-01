#!/bin/bash

#module load cuda

#module load bullxmpi/bullxmpi_1.2.4.1
module load impi/4.1.1.036
#module load impi/4.1.0.030
#module load openmpi/1.8.1-icc


#module load intel/14.0.1
module load intel/15.2.164

module list

cd solver
./compile-anselm-system-intel-mpi-devel.sh

cd ../mesh
./compile-anselm-system-intel-mpi.sh

cd ../permoncube
./compile-anselm-system-intel-mpi.sh

cd ../app
./compile-anselm-system-intel-mpi.sh

