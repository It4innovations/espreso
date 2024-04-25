
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TRIANGLE3", 2, 2, 2, 2, 20, 20, "TOTAL_FETI", "ORTHOGONAL", 24 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for method, max_it in [ (" TOTAL_FETI", 24), ("HYBRID_FETI", 26) ]:
            for projector in [ "             ORTHOGONAL", "ORTHOGONAL_WITH_FACTORS" ]:
                yield run, method, projector, max_it

def run(method, projector, max_it):
    ESPRESOTest.args[7] = method
    ESPRESOTest.args[8] = projector
    ESPRESOTest.args[9] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
