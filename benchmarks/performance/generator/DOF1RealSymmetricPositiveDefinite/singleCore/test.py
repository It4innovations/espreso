
import os, math, unittest, itertools
from estest import ESPRESOTest

results = list()

functions = [
    ("B1 ROWS", "avg"),
    ("K+ SURFACE", "avg"),
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

    ("MESH MEMORY FOOTPRINT [GB]", "size"),
    ("PHYSICAL SOLVER MEMORY FOOTPRINT [GB]", "size"),
    ("FETI SOLVER MEMORY FOOTPRINT [GB]", "size"),
    ("TOTAL MEMORY PER NODE [GB]", "size"),
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
    ESPRESOTest.parse(results, "K+ SURFACE", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ NNZ", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ FACTORS NNZ", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ SOLVER MEMORY [MB]", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "F MEMORY [MB]", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "ITERATIONS TOTAL", ["iterations"])
    ESPRESOTest.parse(results, "MESH MEMORY FOOTPRINT [GB]", ["size"])
    ESPRESOTest.parse(results, "PHYSICAL SOLVER MEMORY FOOTPRINT [GB]", ["size"])
    ESPRESOTest.parse(results, "FETI SOLVER MEMORY FOOTPRINT [GB]", ["size"])
    ESPRESOTest.parse(results, "TOTAL MEMORY PER NODE [GB]", ["size"])
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
        report()
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

def run2D(etype, ex, ey, uniform_domains, method, partial_dual):
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
        for etype, mult in [ ("HEXA8", 1), ("HEXA20", 2), ("TETRA4", 1), ("TETRA10", 2), ("PRISMA6", 1), ("PRISMA15", 2), ("PYRAMID5", 1), ("PYRAMID13", 2) ]:
            ex = ey = ez = 4
            while ex * ey * ez < 66000:
                for uniform_domains, method, partial_dual in testcases:
                    yield run3D, etype, ex // mult, ey // mult, ez // mult, uniform_domains, method, partial_dual
                if ex == ey and ey == ez:
                    ez = 2 * ez
                else:
                    if ex == ey:
                        ey = 2 * ey
                    else:
                        ex = 2 * ex

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