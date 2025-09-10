
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 2, 2, 1, 2, 2, 4, 10, 2, 2, "TOTAL_FETI", "ORTHOGONAL", "DEFAULT", 23 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        yield run, " TOTAL_FETI", "ORTHOGONAL", "DEFAULT", 160
        yield run, " TOTAL_FETI", " CONJUGATE", "DEFAULT",  50
        yield run, "HYBRID_FETI", "ORTHOGONAL", "DEFAULT",  65
        yield run, "HYBRID_FETI", " CONJUGATE", "DEFAULT",  45

def run(method, projector, opt, max_it):
    ESPRESOTest.args[10] = method
    ESPRESOTest.args[11] = projector
    ESPRESOTest.args[12] = opt
    ESPRESOTest.args[13] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")

run.assembler = 1
run.correcness = 1
run.generator = 1
run.feti = 1