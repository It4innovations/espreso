
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 2, 2, 1, 2, 2, 2, 4, 4, 8, "ORTHOGONAL", 96 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for projector, max_it in [ ("ORTHOGONAL     ", 100), ("ORTHOGONAL_FULL", 90) ]:
            yield run, projector, max_it

def run(projector, max_it):
    ESPRESOTest.args[10] = projector
    ESPRESOTest.args[11] = max_it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
