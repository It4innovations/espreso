
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TRIANGLE3", 2, 2, 2, 2, 40, 4, "TOTAL_FETI", "IMPLICIT" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for method in [ "TOTAL_FETI" ]:
            for dualop in [ "    IMPLICIT", "    EXPLICIT", "EXPLICIT_GPU" ]:
                yield run, method, dualop

def run(method, dualop):
    ESPRESOTest.args[7] = method
    ESPRESOTest.args[8] = dualop
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
