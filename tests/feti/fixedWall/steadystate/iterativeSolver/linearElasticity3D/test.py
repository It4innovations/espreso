
import os, unittest
from estest import ESPRESOTest

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 2, 2, 1, 2, 2, 2, 4, 4, 8 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for solver in [ "SMALBE" ]:
            yield run, solver

def run(solver):
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
