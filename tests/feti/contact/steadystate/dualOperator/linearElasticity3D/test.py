
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 2, 2, 1, 2, 2, 2, 4, 4, 8, "TOTAL_FETI" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for method, max_it in [ (" TOTAL_FETI", 93), ("HYBRID_FETI", 96) ]:
            yield run, method, max_it

def run(method, max_it):
    ESPRESOTest.args[10] = method
    ESPRESOTest.args[11] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
