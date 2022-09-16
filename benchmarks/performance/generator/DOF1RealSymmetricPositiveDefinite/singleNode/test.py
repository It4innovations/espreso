
import os, math, unittest, itertools, multiprocessing
from estest import ESPRESOTest

results = list()

functions = [
    ("MPI_COMM_WORLD", "size"),
    ("OMP_NUM_THREADS", "size"),
    ("DOMAINS TOTAL", "size"),
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

    ("K+ SOLVER MEMORY [MB]", "avg"),
    ("F MEMORY [MB]", "avg"),
    ("MESH MEMORY FOOTPRINT [GB]", "size"),
    ("PHYSICAL SOLVER MEMORY FOOTPRINT [GB]", "size"),
    ("FETI SOLVER MEMORY FOOTPRINT [GB]", "size"),
    ("TOTAL MEMORY PER NODE [GB]", "size"),

    ("B1 ROWS", "min"),
    ("B1 ROWS", "max"),
    ("K+ SOLVER MEMORY [MB]", "min"),
    ("K+ SOLVER MEMORY [MB]", "max"),
    ("F MEMORY [MB]", "min"),
    ("F MEMORY [MB]", "max"),
    ("FETI: TFETI SYMBOLIC FACTORIZATION", "min"),
    ("FETI: TFETI SYMBOLIC FACTORIZATION", "max"),
    ("FETI: TFETI NUMERICAL FACTORIZATION", "min"),
    ("FETI: TFETI NUMERICAL FACTORIZATION", "max"),
    ("FETI: TFETI (B * K+ * B') ASSEMBLED", "min"),
    ("FETI: TFETI (B * K+ * B') ASSEMBLED", "max"),
    ("cpg: apply F [~]", "min"),
    ("cpg: apply F [~]", "max"),
    ("cpg TOTAL DURATION", "min"),
    ("cpg TOTAL DURATION", "max"),
    ("FETI[SOLVE] TOTAL DURATION", "min"),
    ("FETI[SOLVE] TOTAL DURATION", "max"),
    ]

testcases = [ (x, y) for x, y in itertools.product([ "UNIFORM", "METIS" ], [ "IMPLICIT_TFETI", "EXPLICIT_TFETI" ]) ]

cores = multiprocessing.cpu_count()
threads = int(os.getenv("OMP_NUM_THREADS"))
mpiprocs = cores // threads

def set(values, size, start, dimension):
    ii = start
    while 1 < size:
        for n in range(2, size + 1):
            if size % n == 0:
                values[ii] = n * values[ii]
                ii = (ii + 1) % dimension
                size = size // n
                break
    return ii

clusters3D = [ 1, 1, 1 ]
domains3D  = [ 1, 1, 1 ]

set(domains3D, threads, set(clusters3D, mpiprocs, 0, 3), 3)
size = max(clusters3D[0] * domains3D[0], max(clusters3D[1] * domains3D[1], clusters3D[2] * domains3D[2]))

domains3D[0] = domains3D[0] * (size // (clusters3D[0] * domains3D[0]))
domains3D[1] = domains3D[1] * (size // (clusters3D[1] * domains3D[1]))
domains3D[2] = domains3D[2] * (size // (clusters3D[2] * domains3D[2]))

clusters2D = [ 1, 1 ]
domains2D  = [ 1, 1 ]

set(domains2D, threads, set(clusters2D, mpiprocs, 0, 2), 2)
size = max(clusters2D[0] * domains2D[0], clusters3D[1] * domains3D[1])

domains2D[0] = domains2D[0] * (size // (clusters2D[0] * domains2D[0]))
domains2D[1] = domains2D[1] * (size // (clusters2D[1] * domains2D[1]))

def report():
    def pick(uniform_domains, method):
        measurement = list()
        for key, values in results:
            if key[2] == uniform_domains and key[3] == method:
                measurement.append(((key[0], key[1]), values))
        return measurement

    ESPRESOTest.parse(results, "MPI_COMM_WORLD", ["size"])
    ESPRESOTest.parse(results, "OMP_NUM_THREADS", ["size"])
    ESPRESOTest.parse(results, "DOMAINS TOTAL", ["size"])
    ESPRESOTest.parse(results, "MESH MEMORY FOOTPRINT [GB]", ["size"])
    ESPRESOTest.parse(results, "PHYSICAL SOLVER MEMORY FOOTPRINT [GB]", ["size"])
    ESPRESOTest.parse(results, "FETI SOLVER MEMORY FOOTPRINT [GB]", ["size"])
    ESPRESOTest.parse(results, "TOTAL MEMORY PER NODE [GB]", ["size"])
    ESPRESOTest.parse(results, "B1 ROWS", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ SURFACE", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ ROWS", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ NNZ", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ FACTORS NNZ", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "K+ SOLVER MEMORY [MB]", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "F MEMORY [MB]", ["avg", "min", "max"])
    ESPRESOTest.parse(results, "ITERATIONS TOTAL", ["iterations"])
    for uniform_domains, method in testcases:
        ESPRESOTest.report_mesurement("_".join([method, uniform_domains, str(mpiprocs), str(threads)]), pick(uniform_domains, method), functions)

class TotalFETIOperator2D(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.processes = mpiprocs
        ESPRESOTest.set_threads(threads)
        ESPRESOTest.ecf = "espreso2D.ecf"
        ESPRESOTest.args = [ "etype", clusters2D[0], clusters2D[1], domains2D[0], domains2D[1], "ex", "ey", "uniform_domains", "TFETI" ]

    def tearDown(self):
        report()
        ESPRESOTest.clean()

    def test_performance(self):
        for etype, mult in [ ("SQUARE4", 1), ("SQUARE8", 2), ("TRIANGLE3", 1), ("TRIANGLE6", 2) ]:
            ex = ey = 8
            while ex * ey < 66000:
                for uniform_domains, method in testcases:
                    yield run2D, etype, ex // mult, ey // mult, uniform_domains, method
                if ex == ey:
                    ey = 2 * ey
                else:
                    ex = 2 * ex

def run2D(etype, ex, ey, uniform_domains, method):
    results.append(((etype, ex * ey, uniform_domains, method), {}))
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[5] = ex - 1
    ESPRESOTest.args[6] = ey - 1
    ESPRESOTest.args[7] = uniform_domains == "UNIFORM"
    ESPRESOTest.args[8] = method
    results[-1] = ((etype, ex * ey, uniform_domains, method), ESPRESOTest.extract(ESPRESOTest.run(), [e[0] for e in functions]))

run2D.performance = 1
run2D.generator = 1

class TotalFETIOperator3D(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.processes = mpiprocs
        ESPRESOTest.set_threads(threads)
        ESPRESOTest.ecf = "espreso3D.ecf"
        ESPRESOTest.args = [ "etype", clusters3D[0], clusters3D[1], clusters3D[2], domains3D[0], domains3D[1], domains3D[2], "ex", "ey", "ez", "uniform_domains", "TFETI" ]

    def tearDown(self):
        report()
        ESPRESOTest.clean()

    def test_performance(self):
        for etype, mult in [ ("HEXA8", 1), ("HEXA20", 2), ("TETRA4", 1), ("TETRA10", 2), ("PRISMA6", 1), ("PRISMA15", 2), ("PYRAMID5", 1), ("PYRAMID13", 2) ]:
            ex = ey = ez = 4
            while ex * ey * ez < 66000:
                for uniform_domains, method in testcases:
                    yield run3D, etype, ex // mult, ey // mult, ez // mult, uniform_domains, method
                if ex == ey and ey == ez:
                    ez = 2 * ez
                else:
                    if ex == ey:
                        ey = 2 * ey
                    else:
                        ex = 2 * ex

def run3D(etype, ex, ey, ez, uniform_domains, method):
    results.append(((etype, ex * ey * ez, uniform_domains, method), {}))
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[7] = ex - 1
    ESPRESOTest.args[8] = ey - 1
    ESPRESOTest.args[9] = ez - 1
    ESPRESOTest.args[10] = uniform_domains == "UNIFORM"
    ESPRESOTest.args[11] = method
    results[-1] = ((etype, ex * ey * ez, uniform_domains, method), ESPRESOTest.extract(ESPRESOTest.run(), [e[0] for e in functions]))

run3D.performance = 1
run3D.generator = 1
