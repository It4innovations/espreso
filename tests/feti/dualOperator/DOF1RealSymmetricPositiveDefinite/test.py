
import os, math, unittest
from estest import ESPRESOTest

testcase = [
    ([ 1, 1, 1,  2, 2, 2,  4, 4, 4 ], 40, 55),
    ([ 1, 1, 1,  4, 4, 4,  4, 4, 4 ], 64, 71),
    ([ 1, 1, 1,  4, 4, 4,  9, 9, 9 ], 71, 87),
    ([ 2, 1, 2,  2, 4, 2,  9, 9, 9 ], 71, 84),
    ([ 2, 2, 2,  2, 2, 2,  9, 9, 9 ], 71, 82),
    ([ 1, 1, 4,  2, 2, 1,  4, 5, 6 ], 86, 99),
    ([ 2, 2, 4,  4, 4, 3,  7, 7, 7 ], 88, 129),
    ]

class TotalFETI(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "cx", "cy", "cz",  "dx", "dy", "dz",  "ex", "ey", "ez", "method", "uniform_clusters" ]
        ESPRESOTest.set_threads(2)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_stability(self):
        for dual in [ "IMPLICIT_TFETI", "EXPLICIT_TFETI" ]:
            for size, uniform, metis in testcase:
                yield run, size, dual, "UNIFORM", uniform
                yield run, size, dual, "METIS", metis

def run(size, dual, uniform_domains, expectation):
    ESPRESOTest.processes = size[0] * size[1] * size[2]
    ESPRESOTest.args[0:9] = size
    ESPRESOTest.args[9] = dual
    ESPRESOTest.args[10] = uniform_domains == "UNIFORM"
    iterations = int(ESPRESOTest.extract(ESPRESOTest.run(), [ "ITERATIONS TOTAL" ])["ITERATIONS TOTAL"])
    if iterations < expectation - 1 or expectation + 1 < iterations:
        ESPRESOTest.raise_error("invalied number of iterations: {0} instead of {1}\n".format(iterations, expectation))

run.stability = 1
run.generator = 1
