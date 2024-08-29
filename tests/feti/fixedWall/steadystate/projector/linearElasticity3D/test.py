
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 2, 2, 1, 2, 2, 2, 4, 4, 8, "ORTHOGONAL", "DEFAULT", 96 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for opt, max_it in [ ("          DEFAULT", 140), ("     WITH_FACTORS", 140), ("WITH_FACTORS|FULL", 110) ]:
            yield run, opt, max_it

def run(opt, max_it):
    ESPRESOTest.args[11] = opt
    ESPRESOTest.args[12] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
