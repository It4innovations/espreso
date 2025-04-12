#!/bin/bash

# assuming Karolina

if [ "$#" -lt 2 ]
then
    echo "not enough arguments"
    exit 1
fi

phase="${1}"
datestr="${2}"

ml mpi4py
ml matplotlib

if [ "${phase}" == "1" ]
then
    mpirun -n 40 python3 benchmarks/dualop_sctria_options/graphs_phase1.py "${datestr}"
    # mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase1.py "${datestr}"
fi
if [ "${phase}" == "2" ]
then
    mpirun -n 8 python3 benchmarks/dualop_sctria_options/graphs_phase2.py "${datestr}"
    # mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase2.py "${datestr}"
fi
if [ "${phase}" == "3" ]
then
    mpirun -n 8 python3 benchmarks/dualop_sctria_options/graphs_phase3.py "${datestr}"
    # mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase3.py "${datestr}"
fi
if [ "${phase}" == "4" ]
then
    mpirun -n 8 python3 benchmarks/dualop_sctria_options/graphs_phase4.py "${datestr}"
    # mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase4.py "${datestr}"
fi
if [ "${phase}" == "5" ]
then
    mpirun -n 8 python3 benchmarks/dualop_sctria_options/graphs_phase5.py "${datestr}"
    # mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase5.py "${datestr}"
fi
if [ "${phase}" == "6" ]
then
    mpirun -n 4 python3 benchmarks/dualop_sctria_options/graphs_phase6.py "${datestr}"
    # mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase6.py "${datestr}"
fi
if [ "${phase}" == "7" ]
then
    mpirun -n 4 python3 benchmarks/dualop_sctria_options/graphs_phase7.py "${datestr}"
    # mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase7.py "${datestr}"
fi
if [ "${phase}" == "8" ]
then
    mpirun -n 4 python3 benchmarks/dualop_sctria_options/graphs_phase8.py "${datestr}"
    # mpirun -n 1 python3 benchmarks/dualop_sctria_options/graphs_phase8.py "${datestr}"
fi
