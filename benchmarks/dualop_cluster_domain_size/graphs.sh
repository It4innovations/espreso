#!/bin/bash

# assuming MN5

ml openmpi/4.1.5-gcc

srun -n 16 -c 1 -G 0 ./benchmarks/dualop_cluster_domain_size/graphs_times_separate.py "$@" &
srun -n 16 -c 1 -G 0 ./benchmarks/dualop_cluster_domain_size/graphs_optima_separate.py "$@" &
srun -n 16 -c 1 -G 0 ./benchmarks/dualop_cluster_domain_size/graphs_times_combined_contour.py "$@" &
srun -n 16 -c 1 -G 0 ./benchmarks/dualop_cluster_domain_size/graphs_times_combined_multiline.py "$@" &
srun -n 16 -c 1 -G 0 ./benchmarks/dualop_cluster_domain_size/graphs_optima_combined.py "$@" &

wait
