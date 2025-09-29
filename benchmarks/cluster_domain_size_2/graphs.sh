#!/bin/bash

# assuming MN5

ml openmpi/4.1.5-gcc

srun -n $((12/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_times_separate_1.py "$@" &
srun -n $((12/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_times_separate_2.py "$@" &
srun -n $((12/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_times_separate_3.py "$@" &
srun -n $((24/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_times_combined_1.py "$@" &
srun -n $((12/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_times_combined_2.py "$@" &
srun -n $((12/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_speedup_1.py "$@" &
srun -n $((12/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_speedup_2.py "$@" &
srun -n $((12/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_speedup_3.py "$@" &
srun -n $((12/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_speedup_4.py "$@" &
srun -n $((12/2)) -c 1 -G 0 ./benchmarks/cluster_domain_size_2/graphs_optima_separate.py "$@" &

wait
