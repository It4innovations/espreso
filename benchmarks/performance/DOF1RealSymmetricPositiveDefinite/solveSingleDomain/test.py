
import os, math, unittest
from estest import ESPRESOTest

results = list()

functions = [
    ("FETI: TFETI SYMBOLIC FACTORIZATION", "avg"),
    ("FETI[SET] TOTAL DURATION", "avg"),
    ("K+ ROWS", "avg"),
    ("K+ NNZ", "avg"),
    ("K+ FACTORS NNZ", "avg"),
    ("FETI: TFETI NUMERICAL FACTORIZATION", "avg"),
    ("FETI: TFETI (B * K+ * B') ASSEMBLED", "avg"),
    ("FETI[UPDATE] TOTAL DURATION", "avg"),
    ("cpg: apply F [~]", "avg"),
    ("cpg TOTAL DURATION", "avg")
    ]

def report():
    def pick(uniform_domains, method, partial_dual):
        measurement = list()
        for key, values in results:
            if key[2] == uniform_domains and key[3] == method and key[4] == partial_dual:
                measurement.append(((key[0], key[1]), values))
        return measurement

    ESPRESOTest.parse(results, "K+ ROWS", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ NNZ", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ FACTORS NNZ", ["avg", "min", "max"])
    ESPRESOTest.report_mesurement("IMPLICIT_UNIFORM_FULL", pick("UNIFORM", "IMPLICIT_TFETI", "FULL"), functions)
    ESPRESOTest.report_mesurement("EXPLICIT_UNIFORM_FULL", pick("UNIFORM", "EXPLICIT_TFETI", "FULL"), functions)
    ESPRESOTest.report_mesurement("EXPLICIT_UNIFORM_PARTIAL", pick("UNIFORM", "EXPLICIT_TFETI", "PARTIAL"), functions)
    ESPRESOTest.report_mesurement("IMPLICIT_METIS_FULL", pick("METIS", "IMPLICIT_TFETI", "FULL"), functions)
    ESPRESOTest.report_mesurement("EXPLICIT_METIS_FULL", pick("METIS", "EXPLICIT_TFETI", "FULL"), functions)
    ESPRESOTest.report_mesurement("EXPLICIT_METIS_PARTIAL", pick("METIS", "EXPLICIT_TFETI", "PARTIAL"), functions)

class TotalFETIOperator2D(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.processes = 1
        ESPRESOTest.ecf = "espreso2D.ecf"
        ESPRESOTest.args = [ "etype", 2, 2, "uniform_domains", "TFETI", "restricted_dual" ]

    def tearDown(self):
        report()
        ESPRESOTest.clean()

    def test_performance(self):
        for etype in [ "SQUARE4", "SQUARE8", "TRIANGLE3", "TRIANGLE6" ]:
            for size in range(7, 14):
                yield run2D, etype, int(math.sqrt(2 ** size) - 1), "UNIFORM", "IMPLICIT_TFETI", "FULL"
                yield run2D, etype, int(math.sqrt(2 ** size) - 1), "UNIFORM", "EXPLICIT_TFETI", "FULL"
                yield run2D, etype, int(math.sqrt(2 ** size) - 1), "UNIFORM", "EXPLICIT_TFETI", "PARTIAL"
                yield run2D, etype, int(math.sqrt(2 ** size) - 1), "METIS", "IMPLICIT_TFETI", "FULL"
                yield run2D, etype, int(math.sqrt(2 ** size) - 1), "METIS", "EXPLICIT_TFETI", "FULL"
                yield run2D, etype, int(math.sqrt(2 ** size) - 1), "METIS", "EXPLICIT_TFETI", "PARTIAL"

class TotalFETIOperator3D(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.processes = 1
        ESPRESOTest.ecf = "espreso3D.ecf"
        ESPRESOTest.args = [ "etype", 2, 2, 2, "uniform_domains", "TFETI", "restricted_dual" ]

    def tearDown(self):
        report()
        ESPRESOTest.clean()

    def test_performance(self):
        for etype in [ "HEXA8", "HEXA20", "TETRA4", "TETRA10", "PRISMA6", "PRISMA15", "PYRAMID5", "PYRAMID13" ]:
            for size in range(7, 14):
                yield run3D, etype, int(math.pow(2 ** size, 1./3) - 1), "UNIFORM", "IMPLICIT_TFETI", "FULL"
                yield run3D, etype, int(math.pow(2 ** size, 1./3) - 1), "UNIFORM", "EXPLICIT_TFETI", "FULL"
                yield run3D, etype, int(math.pow(2 ** size, 1./3) - 1), "UNIFORM", "EXPLICIT_TFETI", "PARTIAL"
                yield run3D, etype, int(math.pow(2 ** size, 1./3) - 1), "METIS", "IMPLICIT_TFETI", "FULL"
                yield run3D, etype, int(math.pow(2 ** size, 1./3) - 1), "METIS", "EXPLICIT_TFETI", "FULL"
                yield run3D, etype, int(math.pow(2 ** size, 1./3) - 1), "METIS", "EXPLICIT_TFETI", "PARTIAL"

def run2D(etype, size, uniform_domains, method, partial_dual):
    ESPRESOTest.set_threads(1)
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[1] = size
    ESPRESOTest.args[2] = size
    ESPRESOTest.args[3] = uniform_domains == "UNIFORM"
    ESPRESOTest.args[4] = method
    ESPRESOTest.args[5] = partial_dual == "PARTIAL"
    results.append(((etype, (size + 1) * (size + 1), uniform_domains, method, partial_dual), ESPRESOTest.extract(ESPRESOTest.run(), [e[0] for e in functions])))

def run3D(etype, size, uniform_domains, method, partial_dual):
    ESPRESOTest.set_threads(1)
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[1] = size
    ESPRESOTest.args[2] = size
    ESPRESOTest.args[3] = size
    ESPRESOTest.args[4] = uniform_domains == "UNIFORM"
    ESPRESOTest.args[5] = method
    ESPRESOTest.args[6] = partial_dual == "PARTIAL"
    results.append(((etype, (size + 1) * (size + 1) * (size + 1), uniform_domains, method, partial_dual), ESPRESOTest.extract(ESPRESOTest.run(), [e[0] for e in functions])))

run2D.performance = 1
run2D.generator = 1
run3D.performance = 1
run3D.generator = 1
