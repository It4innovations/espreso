#!/usr/bin/python

import os
import sys

os.environ["LD_LIBRARY_PATH"] += ":./libs"

if "mesh" in sys.argv:
    if "valgrind" in sys.argv:
        os.system("valgrind ./build/mesh/src/esmesh")
    else:
        os.system("./build/mesh/src/esmesh")
else:
    os.system("./espreso")



