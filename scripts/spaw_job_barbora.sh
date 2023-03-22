#!/bin/bash

qsub -q qcpu -N espreso_barbora_cpu -l select=1:ncpus=36,walltime=48:00:00,msr=1.4.0 -A OPEN-22-7 ./SC_paper_data_barbora_cpu.sh
qsub -q qgpu -N espreso_barbora_gpu -l select=1:ncpus=24,walltime=48:00:00,msr=1.4.0 -A OPEN-22-7 ./SC_paper_data_barbora_gpu.sh