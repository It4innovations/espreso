
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TRIANGLE3", 2, 2, 2, 2, 20, 20, "TOTAL_FETI", "ORTHOGONAL", "DEFAULT", 24 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for opt in [ "          DEFAULT", "     WITH_FACTORS", "             FULL", "WITH_FACTORS|FULL" ]:
            yield run, " TOTAL_FETI", "ORTHOGONAL", opt, 24
        for opt in [ "          DEFAULT", "     WITH_FACTORS"]:
            yield run, "HYBRID_FETI", "ORTHOGONAL", opt, 26
        for opt in [ "             FULL", "WITH_FACTORS|FULL" ]:
            yield run, "HYBRID_FETI", "ORTHOGONAL", opt, 22

def run(method, projector, opt, max_it):
    ESPRESOTest.args[ 7] = method
    ESPRESOTest.args[ 8] = projector
    ESPRESOTest.args[ 9] = opt
    ESPRESOTest.args[10] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
