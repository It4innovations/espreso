
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TETRA4", 2, 2, 1, 2, 2, 2, 4, 4, 8, "TOTAL_FETI", "ORTHOGONAL", 23 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for method, max_it in [ (" TOTAL_FETI", 23), ("HYBRID_FETI", 25) ]:
            for projector in [ "             ORTHOGONAL", "ORTHOGONAL_WITH_FACTORS" ]:
                yield run, method, projector, max_it

def run(method, projector, max_it):
    ESPRESOTest.args[10] = method
    ESPRESOTest.args[11] = projector
    ESPRESOTest.args[12] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
