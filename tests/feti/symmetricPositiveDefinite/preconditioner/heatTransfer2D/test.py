
import os, unittest
from estest import ESPRESOTest

cases = [
    ("     NONE", 85),
    ("   LUMPED", 55),
    ("DIRICHLET", 23),
    ]

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TRIANGLE3", 2, 2, 2, 2, 20, 20, "NONE", "50" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for prec, iterations in cases:
            yield run, prec, iterations

def run(prec, iterations):
    ESPRESOTest.args[7] = prec
    ESPRESOTest.args[8] = iterations
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
