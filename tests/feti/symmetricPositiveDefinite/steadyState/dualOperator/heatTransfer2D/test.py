
import os, unittest
from estest import ESPRESOTest

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TRIANGLE3", 2, 2, 2, 2, 20, 20, "TOTAL_FETI", "IMPLICIT" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for dualop in [ "    IMPLICIT", "    EXPLICIT", "EXPLICIT_GPU" ]:
            yield run, " TOTAL_FETI", dualop
        for dualop in [ "    IMPLICIT" ]:
            yield run, "HYBRID_FETI", dualop

def run(method, dualop):
    ESPRESOTest.args[7] = method
    ESPRESOTest.args[8] = dualop
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
