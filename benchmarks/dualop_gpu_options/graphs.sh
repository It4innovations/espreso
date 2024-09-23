#/bin/bash

# assuming Karolina

ml mpi4py
ml matplotlib

mpirun -n 8 python3 benchmarks/dualop_gpu_options/graphs.py
