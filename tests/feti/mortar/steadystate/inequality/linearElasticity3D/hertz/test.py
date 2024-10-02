
import os, unittest
from estest import ESPRESOTest

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 2, 5, "TETRA4", 2, 9, 230 ]
        ESPRESOTest.processes = 2
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        iterations = [ 230, 290, 340, 360 ]
        it = 0
        for el1 in [ " HEXA8", "TETRA4" ]:
            for el2 in [ " HEXA8", "TETRA4" ]:
                yield run, el1, el2, iterations[it]
                it = it + 1

def run(el1, el2, it):
    ESPRESOTest.args[0] = el1
    ESPRESOTest.args[3] = el2
    ESPRESOTest.args[6] = it
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("{0}_{1}.emr".format(el1, el2))
