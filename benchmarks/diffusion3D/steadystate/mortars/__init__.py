
from estest import ESPRESOTest

def setup():
    ESPRESOTest.processes = 2
    ESPRESOTest.env["OMP_NUM_THREADS"] = "4"
    ESPRESOTest.env["SOLVER_NUM_THREADS"] = "4"
    ESPRESOTest.env["PAR_NUM_THREADS"] = "4"