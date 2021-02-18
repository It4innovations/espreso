#!/bin/bash

ml Forge

. env/modules.barbora.icpc

#ml Python/3.8.2-GCCcore-9.3.0

./waf configure -m debug  --solver=cuda

./waf

. env/threading.default 1

mpirun -n 4 espreso 8 2 2 2 2 3 3  -vvv 
 
#ddt -n 1 espreso 8 1 1 1 1 3 3  -vvv 
