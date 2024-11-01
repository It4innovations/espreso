
import os, unittest
from estest import ESPRESOTest

cases = [
    ("          PCG", "     NONE", 85),
    ("          PCG", "DIRICHLET", 27),
    ("       SMALBE", "     NONE", 85),
    ("       SMALBE", "DIRICHLET", 27),
    ("orthogonalPCG", "     NONE", 63),
    ("orthogonalPCG", "DIRICHLET", 25),
    ]

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TRIANGLE3", 2, 2, 2, 2, 20, 20, "PCG", "NONE", "50" ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for solver, prec, iterations in cases:
            yield run, solver, prec, iterations

def run(solver, prec, iterations):
    ESPRESOTest.args[7] = solver
    ESPRESOTest.args[8] = prec
    ESPRESOTest.args[9] = iterations
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
