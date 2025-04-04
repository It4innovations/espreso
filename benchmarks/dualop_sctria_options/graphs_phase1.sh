#!/bin/bash

# assuming Karolina

ml mpi4py
ml matplotlib

mpirun -n 40 python3 benchmarks/dualop_sctria_options/graphs_phase1.py "$@"
# mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase1.py "$@"
