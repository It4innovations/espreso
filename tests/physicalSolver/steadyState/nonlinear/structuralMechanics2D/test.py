
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "SQUARE4", 2, 2, 2, 2, 5, 50, "MKLPDSS" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_physical_solver(self):
        for solver in [ "MKLPDSS" ]:
            yield run, solver

def run(solver):
    ESPRESOTest.args[7] = solver
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")

run.physical_solver = 1
run.generator = 1