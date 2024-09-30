
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "etype", 2, 2, 1, 4, 4, 2, 8, 2, 6 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_correctness(self):
        for etype in [ "    HEXA8", "   TETRA4", "  PRISMA6", " PYRAMID5" ]:
        # for etype in [ "    HEXA8", "   HEXA20", "   TETRA4", "  TETRA10", "  PRISMA6", " PRISMA15", " PYRAMID5", "PYRAMID13" ]:
            yield run, etype

def run(etype):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))

run.assembler = 1
run.correcness = 1
run.generator = 1
