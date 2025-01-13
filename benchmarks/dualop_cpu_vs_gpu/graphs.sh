#!/bin/bash

# assuming Karolina

ml mpi4py
ml matplotlib

summarize_datestr=20250113_233607
datestr="$(date +%Y%m%d_%H%M%S)"

mpirun -n 6 python3 benchmarks/dualop_cpu_vs_gpu/graphs_amort.py "${summarize_datestr}" "${datestr}" &
mpirun -n 48 python3 benchmarks/dualop_cpu_vs_gpu/graphs_time_machines.py "${summarize_datestr}" "${datestr}" &
mpirun -n 16 python3 benchmarks/dualop_cpu_vs_gpu/graphs_time_dualops.py "${summarize_datestr}" "${datestr}" &
mpirun -n 8 python3 benchmarks/dualop_cpu_vs_gpu/graphs_time_all.py "${summarize_datestr}" "${datestr}" &
mpirun -n 6 python3 benchmarks/dualop_cpu_vs_gpu/graphs_bestdualop.py "${summarize_datestr}" "${datestr}" &
mpirun -n 6 python3 benchmarks/dualop_cpu_vs_gpu/graphs_speedup_xndofs.py "${summarize_datestr}" "${datestr}" &
mpirun -n 6 python3 benchmarks/dualop_cpu_vs_gpu/graphs_speedup_xniters.py "${summarize_datestr}" "${datestr}" &
mpirun -n 48 python3 benchmarks/dualop_cpu_vs_gpu/graphs_time_ndofs.py "${summarize_datestr}" "${datestr}" &

wait
