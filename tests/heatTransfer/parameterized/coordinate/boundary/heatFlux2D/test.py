
import os, unittest
from estest import ESPRESOTest

class Assembler(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "etype", 2, 2, 2, 2, 20, 20 ]
        ESPRESOTest.processes = 4
        ESPRESOTest.set_threads(2)

    def teardown():
        ESPRESOTest.clean()

    def test_correctness(self):
        for etype in [ "SQUARE4", "SQUARE8", "TRIANGLE3", "TRIANGLE6" ]:
            yield run, etype

def run(etype):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))

run.assembler = 1
run.correcness = 1
run.generator = 1