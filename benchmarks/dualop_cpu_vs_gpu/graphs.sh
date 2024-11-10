#/bin/bash

# assuming Karolina

ml mpi4py
ml matplotlib

mpirun -n 19 python3 benchmarks/dualop_cpu_vs_gpu/graphs.py
