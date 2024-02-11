
import os, unittest
from estest import ESPRESOTest

cases = [
    ("MPRGP", "     NONE", 40),
    ]

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 2, 2, 1, 2, 2, 2, 4, 4, 8, "MPRGP", "NONE", 980 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for solver, prec, iterations in cases:
            yield run, solver, prec, iterations

def run(solver, prec, iterations):
    ESPRESOTest.args[10] = solver
    ESPRESOTest.args[11] = prec
    ESPRESOTest.args[12] = iterations
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
