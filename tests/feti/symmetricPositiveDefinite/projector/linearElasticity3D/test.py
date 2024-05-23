
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TETRA4", 2, 2, 1, 2, 2, 2, 4, 4, 8, "TOTAL_FETI", "ORTHOGONAL", "DEFAULT", 15 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for opt in [ "          DEFAULT", "     WITH_FACTORS", "             FULL", "WITH_FACTORS|FULL" ]:
            yield run, " TOTAL_FETI", "ORTHOGONAL", opt, 15
        for opt in [ "          DEFAULT", "     WITH_FACTORS"]:
            yield run, "HYBRID_FETI", "ORTHOGONAL", opt, 17
        for opt in [ "             FULL", "WITH_FACTORS|FULL" ]:
            yield run, "HYBRID_FETI", "ORTHOGONAL", opt, 13

def run(method, projector, opt, max_it):
    ESPRESOTest.args[10] = method
    ESPRESOTest.args[11] = projector
    ESPRESOTest.args[12] = opt
    ESPRESOTest.args[13] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
