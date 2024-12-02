
import os, unittest
from estest import ESPRESOTest

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TETRA4", 2, 2, 1, 2, 2, 2, 4, 4, 8, "TOTAL_FETI", "IMPLICIT", 20 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for dualop in [ "    IMPLICIT", "    EXPLICIT", "EXPLICIT_GPU" ]:
            yield run, " TOTAL_FETI", dualop, 15
        for dualop in [ "    IMPLICIT" ]:
            yield run, "HYBRID_FETI", dualop, 20

def run(method, dualop, max_it):
    ESPRESOTest.args[10] = method
    ESPRESOTest.args[11] = dualop
    ESPRESOTest.args[12] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
