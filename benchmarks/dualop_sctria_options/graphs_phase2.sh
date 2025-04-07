#!/bin/bash

# assuming Karolina

ml mpi4py
ml matplotlib

mpirun -n 8 python3 benchmarks/dualop_sctria_options/graphs_phase2.py "$@"
# mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase2.py "$@"
