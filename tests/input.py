
from utils import *
import unittest
import glob
import string

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(ESPRESO_TESTS)

class ESPRESOInput(unittest.TestCase):

    espreso = Espreso()

    def openfoam(self, path, file):
        config = { "ENV::TESTING_LEVEL": 0, "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULTS": 0 }
        self.espreso.run(len(glob.glob(path + "/processor*")), path, config)

    def workbench(self, path, file):
        config = { "ENV::TESTING_LEVEL": 0, "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULTS": 0 }
        self.espreso.run(1, path, config)

if __name__ == '__main__':

    openfoam  = os.path.join(ESPRESO_TESTS, "examples", "input", "openfoam")
    workbench = os.path.join(ESPRESO_TESTS, "examples", "input", "workbench")

    for name, path, file in TestCaseCreator.gather(openfoam, ".ecf"):
        TestCaseCreator.create_test(ESPRESOInput, ESPRESOInput.openfoam, "openfoam_" + name, path, file)

    for name, path, file in TestCaseCreator.gather(workbench, ".ecf"):
        TestCaseCreator.create_test(ESPRESOInput, ESPRESOInput.workbench, "workbench_" + name, path, file)

    unittest.main()
