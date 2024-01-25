
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "SQUARE4", 2, 2, 2, 2, 20, 20, "FETI" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for solver in [ "MKLPDSS", "   FETI" ]:
            yield run, solver

def run(solver):
    ESPRESOTest.args[7] = solver
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
