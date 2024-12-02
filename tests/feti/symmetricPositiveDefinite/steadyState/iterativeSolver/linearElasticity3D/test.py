
import os, unittest
from estest import ESPRESOTest

cases = [
    ("          PCG", "     NONE", 700),
    ("          PCG", "DIRICHLET", 190),
    ("       SMALBE", "     NONE", 800), # there are slightly different convergence criteria 
    ("       SMALBE", "DIRICHLET", 190),
    ("orthogonalPCG", "     NONE", 290),
    ("orthogonalPCG", "DIRICHLET", 120),
    ]

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "TETRA4", 2, 2, 1, 2, 2, 2, 4, 4, 8, "PCG", "NONE", 200 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        for solver, prec, iterations in cases:
            yield run, solver, prec, iterations

def run(solver, prec, iterations):
    ESPRESOTest.args[10] = solver
    ESPRESOTest.args[11] = prec
    ESPRESOTest.args[12] = iterations
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
