
import os, unittest
from estest import ESPRESOTest

cases = [
    ("     NONE", 700),
    ("   LUMPED", 290),
    ("DIRICHLET", 190),
    ]

class Assembler(unittest.TestCase):

    def startTest(self, event):
        event.test._testFunc = dummy
        event.handled = True

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TETRA4", 2, 2, 1, 2, 2, 2, 4, 4, 8, "NONE", 200 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for prec, iterations in cases:
            yield run, prec, iterations

def run(prec, iterations):
    ESPRESOTest.args[10] = prec
    ESPRESOTest.args[11] = iterations
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
