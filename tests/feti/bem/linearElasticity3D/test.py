
import os, unittest
from estest import ESPRESOTest

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ 2, 2, 1, 2, 2, 2, 3, 3, 5, "FEM_LOADED" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for method in [ "FEM_LOADED", "BEM       " ]:
            yield run, method

def run(method):
    ESPRESOTest.args[9] = method
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")