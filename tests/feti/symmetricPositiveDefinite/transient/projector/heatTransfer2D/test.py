
import os, unittest
from estest import ESPRESOTest

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TRIANGLE3", 2, 2, 8, 8, 5, 5, "TOTAL_FETI", "ORTHOGONAL", "DEFAULT", 24 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        yield run, " TOTAL_FETI", "ORTHOGONAL", "DEFAULT", 45
        yield run, " TOTAL_FETI", " CONJUGATE", "DEFAULT", 25

def run(method, projector, opt, max_it):
    ESPRESOTest.args[ 7] = method
    ESPRESOTest.args[ 8] = projector
    ESPRESOTest.args[ 9] = opt
    ESPRESOTest.args[10] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
