
from estest import ESPRESOTest

def setup():
    ESPRESOTest.processes = 1
    ESPRESOTest.env["OMP_NUM_THREADS"] = "1"
    ESPRESOTest.env["SOLVER_NUM_THREADS"] = "1"
    ESPRESOTest.env["PAR_NUM_THREADS"] = "1"