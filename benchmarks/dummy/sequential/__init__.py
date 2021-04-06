
from estest import ESPRESOTest

def setup():
    ESPRESOTest.mpirun = [ "mpirun", "-n" ]
    ESPRESOTest.processes = 1
    ESPRESOTest.env["OMP_NUM_THREADS"] = "1"
    ESPRESOTest.env["SOLVER_NUM_THREADS"] = "1"
    ESPRESOTest.env["PAR_NUM_THREADS"] = "1"

def teardown():
    ESPRESOTest.mpirun = [ "mpirun", "--map-by", "socket:pe=2", "--bind-to", "core", "-n" ]
