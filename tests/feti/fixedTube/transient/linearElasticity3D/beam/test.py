
import os, unittest
from estest import ESPRESOTest

class Solver(unittest.TestCase):

    def setUp(self):
        ESPRESOTest.path = os.path.dirname(__file__)
        ESPRESOTest.args = [ "HEXA8", 1, 1, 1, 15, 1, 1, 10, 5, 4 ]
        ESPRESOTest.processes = 1
        ESPRESOTest.set_threads(8)

    def tearDown(self):
        ESPRESOTest.clean()

    def test_feti(self):
        ESPRESOTest.run()
        ESPRESOTest.compare_emr("espreso.emr")
