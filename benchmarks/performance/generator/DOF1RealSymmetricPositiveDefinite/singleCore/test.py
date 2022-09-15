
import os, math, unittest, itertools
from estest import ESPRESOTest

results = list()

functions = [
    ("B1 ROWS", "avg"),
    ("K+ ROWS", "avg"),
    ("K+ NNZ", "avg"),
    ("K+ FACTORS NNZ", "avg"),
    ("FETI: TFETI SYMBOLIC FACTORIZATION", "avg"),
    ("FETI: TFETI NUMERICAL FACTORIZATION", "avg"),
    ("FETI: TFETI (B * K+ * B') ASSEMBLED", "avg"),
    ("cpg: apply F [~]", "avg"),
    ("ITERATIONS TOTAL", "iterations"),
    ("cpg TOTAL DURATION", "avg"),
    ("FETI[SET] TOTAL DURATION", "avg"),
    ("FETI[UPDATE] TOTAL DURATION", "avg"),
    ("FETI[SOLVE] TOTAL DURATION", "avg"),
    ]

testcases = [ (x, y, z) for x, y, z in itertools.product([ "UNIFORM", "METIS" ], [ "IMPLICIT_TFETI", "EXPLICIT_TFETI" ], [ "FULL", "PARTIAL" ]) ]

def report():
    def pick(uniform_domains, method, partial_dual):
        measurement = list()
        for key, values in results:
            if key[2] == uniform_domains and key[3] == method and key[4] == partial_dual:
                measurement.append(((key[0], key[1]), values))
        return measurement

    ESPRESOTest.parse(results, "B1 ROWS", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ ROWS", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ NNZ", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ FACTORS NNZ", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "ITERATIONS TOTAL", ["iterations"])
    for uniform_domains, method, partial_dual in testcases:
        ESPRESOTest.report_mesurement("_".join([method, partial_dual, uniform_domains]), pick(uniform_domains, method, partial_dual), functions)

class TotalFETIOperator2D(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.processes = 1
        ESPRESOTest.set_threads(1)
        ESPRESOTest.ecf = "espreso2D.ecf"
        ESPRESOTest.args = [ "etype", "ex", "ey", "uniform_domains", "TFETI", "restricted_dual" ]

    def tearDown(self):
        # report()
        ESPRESOTest.clean()

    def test_performance(self):
        for etype, mult in [ ("SQUARE4", 1), ("SQUARE8", 2), ("TRIANGLE3", 1), ("TRIANGLE6", 2) ]:
            ex = ey = 8
            while ex * ey < 66000:
                for uniform_domains, method, partial_dual in testcases:
                    yield run2D, etype, ex // mult, ey // mult, uniform_domains, method, partial_dual
                if ex == ey:
                    ey = 2 * ey
                else:
                    ex = 2 * ex
                print(ex, ey)

def run2D(etype, ex, ey, uniform_domains, method, partial_dual):
    return
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[1] = ex - 1
    ESPRESOTest.args[2] = ey - 1
    ESPRESOTest.args[3] = uniform_domains == "UNIFORM"
    ESPRESOTest.args[4] = method
    ESPRESOTest.args[5] = partial_dual == "PARTIAL"
    results.append(((etype, ex * ey, uniform_domains, method, partial_dual), ESPRESOTest.extract(ESPRESOTest.run(), [e[0] for e in functions])))

run2D.performance = 1
run2D.generator = 1

class TotalFETIOperator3D(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.processes = 1
        ESPRESOTest.set_threads(1)
        ESPRESOTest.ecf = "espreso3D.ecf"
        ESPRESOTest.args = [ "etype", "ex", "ey", "ez", "uniform_domains", "TFETI", "restricted_dual" ]

    def tearDown(self):
        report()
        ESPRESOTest.clean()

    def test_performance(self):
        # for etype in [ "HEXA8", "HEXA20", "TETRA4", "TETRA10", "PRISMA6", "PRISMA15", "PYRAMID5", "PYRAMID13" ]:
        for etype in [ "HEXA8" ]:
            ex = ey = ez = 4
            while ex * ey * ez < 66000:
                print(testcases)
                for uniform_domains, method, partial_dual in testcases:
                    yield run3D, etype, ex, ey, ez, uniform_domains, method, partial_dual
                if ex == ey and ey == ez:
                    ez = 2 * ez
                else:
                    if ex == ey:
                        ey = 2 * ey
                    else:
                        ex = 2 * ex
                print(ex, ey, ez)

def run3D(etype, ex, ey, ez, uniform_domains, method, partial_dual):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[1] = ex - 1
    ESPRESOTest.args[2] = ey - 1
    ESPRESOTest.args[3] = ez - 1
    ESPRESOTest.args[4] = uniform_domains == "UNIFORM"
    ESPRESOTest.args[5] = method
    ESPRESOTest.args[6] = partial_dual == "PARTIAL"
    results.append(((etype, ex * ey * ez, uniform_domains, method, partial_dual), ESPRESOTest.extract(ESPRESOTest.run(), [e[0] for e in functions])))

run3D.performance = 1
run3D.generator = 1
