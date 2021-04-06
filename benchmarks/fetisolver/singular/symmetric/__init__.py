
from estest import ESPRESOTest

def setup():
    ESPRESOTest.processes = 6
    ESPRESOTest.env["OMP_NUM_THREADS"] = "2"
    ESPRESOTest.env["SOLVER_NUM_THREADS"] = "2"
    ESPRESOTest.env["PAR_NUM_THREADS"] = "2"