
import os, unittest
from estest import ESPRESOTest

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 2, 2, 1, 1, 1, 2, 4, 4, 12 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for etype in [ "    HEXA8", "   HEXA20", "   TETRA4", "  TETRA10", "  PRISMA6", " PRISMA15", " PYRAMID5", "PYRAMID13" ]:
            yield run, etype

def run(etype):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("{}.emr".format(etype))

run.feti = 1