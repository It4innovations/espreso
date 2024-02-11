
import os, unittest
from estest import ESPRESOTest

class PhysicalSolver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 2, 2, 1, 2, 2, 2, 4, 4, 8, "FETI" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_physical_solver(self):
        for solver in [ "   FETI" ]:
            yield run, solver

def run(solver):
    ESPRESOTest.args[10] = solver
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")

run.physical_solver = 1
run.generator = 1