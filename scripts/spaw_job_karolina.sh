#!/bin/bash
 
qsub -q qcpu -N espreso_karolina_cpu -l select=1:ncpus=128,walltime=48:00:00,msr=1.5.0 -A OPEN-22-7 SC_paper_data.sh karolina cpu
qsub -q qgpu -N espreso_karolina_gpu -l select=1:ncpus=128,walltime=48:00:00,msr=1.5.0 -A OPEN-22-7 SC_paper_data.sh karolina gpu