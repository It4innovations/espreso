
import os, math, unittest, multiprocessing, glob
from estest import ESPRESOTest

tab = list()
functions = [
    ("FETI: TFETI SYMBOLIC FACTORIZATION", "min"),
    ("FETI: TFETI SYMBOLIC FACTORIZATION", "avg"),
    ("FETI: TFETI SYMBOLIC FACTORIZATION", "max"),
    ("FETI[SET] TOTAL DURATION", "avg"),
    ("FETI: TFETI NUMERICAL FACTORIZATION", "min"),
    ("FETI: TFETI NUMERICAL FACTORIZATION", "avg"),
    ("FETI: TFETI NUMERICAL FACTORIZATION", "max"),
    ("FETI[UPDATE] TOTAL DURATION", "avg"),
    ("cpg: apply F [~]", "min"),
    ("cpg: apply F [~]", "avg"),
    ("cpg: apply F [~]", "max"),
    ("cpg TOTAL DURATION", "avg")
    ]

class SingleNodePerformance(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "etype", 1, 1, 2, 2, 5, 5 ]
    
    def tearDown(self):
        ESPRESOTest.report_mesurement(tab, functions)
        ESPRESOTest.clean()

    def test_performance(self):
        for filename in glob.glob(os.path.dirname(__file__) + "/run.*.sh"):
            prefix, system, mpilib, suffix = os.path.split(filename)[1].split(".")
            if mpilib == ESPRESOTest.mpilib and system == ESPRESOTest.system:
                file = open(filename, "r")
                for line in file:
                    mpirun, procs, args = ESPRESOTest.parse_run_command(line)
                    for etype in [ "SQUARE4" ]: #, "SQUARE8", "TRIANGLE3", "TRIANGLE6" ]:
                        for size in range(7, 8):
                            e = math.sqrt(2 ** size) - 1
                            yield run, int(procs), etype, int(args[1]), int(args[2]), int(args[3]), int(args[4]), int(e), int(e), mpirun

def run(procs, etype, cx, cy, dx, dy, ex, ey, mpirun):
    ESPRESOTest.mpirun = mpirun
    ESPRESOTest.processes = procs
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[1] = cx
    ESPRESOTest.args[2] = cy
    ESPRESOTest.args[3] = dx
    ESPRESOTest.args[4] = dy
    ESPRESOTest.args[5] = ex
    ESPRESOTest.args[6] = ey
    tab.append(((etype, procs, (ex + 1) * (ex + 1)), ESPRESOTest.extract(ESPRESOTest.run(), [e[0] for e in functions])))

run.performance = 1
run.generator = 1
