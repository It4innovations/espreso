
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TRIANGLE3", 2, 2, 2, 2, 40, 4, "TOTAL_FETI", "IMPLICIT", 80 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for dualop in [ "    IMPLICIT", "    EXPLICIT", "EXPLICIT_GPU" ]:
            yield run, " TOTAL_FETI", dualop, 75
        for dualop in [ "    IMPLICIT" ]:
            yield run, "HYBRID_FETI", dualop, 80

def run(method, dualop, max_it):
    ESPRESOTest.args[7] = method
    ESPRESOTest.args[8] = dualop
    ESPRESOTest.args[9] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
