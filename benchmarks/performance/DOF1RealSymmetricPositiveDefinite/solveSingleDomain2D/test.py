
import os, math, unittest
from estest import ESPRESOTest

tab = list()
functions = [
    ("FETI: TFETI SYMBOLIC FACTORIZATION", "avg"),
    ("FETI[SET] TOTAL DURATION", "avg"),
    ("FETI: TFETI NUMERICAL FACTORIZATION", "avg"),
    ("FETI[UPDATE] TOTAL DURATION", "avg"),
    ("cpg: apply F [~]", "avg"),
    ("cpg TOTAL DURATION", "avg")
    ]

class SequentialPerformance(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.processes = 1
        ESPRESOTest.args = [ "etype", 2, 2 ]
    
    def tearDown(self):
        ESPRESOTest.report_mesurement(tab, functions)
        ESPRESOTest.clean()

    def test_performance(self):
        for etype in [ "SQUARE4", "SQUARE8", "TRIANGLE3", "TRIANGLE6" ]:
            for size in range(7, 15):
                yield run, etype, int(math.sqrt(2 ** size) - 1)

def run(etype, size):
    ESPRESOTest.set_threads(1)
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[1] = size
    ESPRESOTest.args[2] = size
    tab.append(((etype, (size + 1) * (size + 1)), ESPRESOTest.extract(ESPRESOTest.run(), [e[0] for e in functions])))

run.performance = 1
run.generator = 1
