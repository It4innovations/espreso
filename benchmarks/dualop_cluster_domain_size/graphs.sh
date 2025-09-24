#!/bin/bash

# assuming MN5

ml openmpi/4.1.5-gcc

srun -n 48 -c 1 -G 0 ./benchmarks/dualop_cluster_domain_size/graphs.py "$@"
# srun -n 1 -c 1 -G 0 ./benchmarks/dualop_cluster_domain_size/graphs.py "$@"
