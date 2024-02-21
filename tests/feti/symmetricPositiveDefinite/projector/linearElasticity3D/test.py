
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TETRA4", 2, 2, 1, 2, 2, 2, 4, 4, 8, "ORTHOGONAL" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for projector in [ "             ORTHOGONAL", "ORTHOGONAL_WITH_FACTORS" ]:
            yield run, projector

def run(projector):
    ESPRESOTest.args[10] = projector
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
