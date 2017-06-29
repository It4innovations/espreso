
import os
import sys

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(ESPRESO_TESTS, "utils"))

from testing import *
import unittest
import glob
import shutil

class ESPRESOInput(unittest.TestCase):

    espreso = Espreso()

    def openfoam(self, path, file):
        config = { "ENV::TESTING_LEVEL": 0, "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULTS": 0 }
        self.espreso.run(len(glob.glob(path + "/processor*")), path, config)

    def workbench(self, path, file):
        config = { "ENV::TESTING_LEVEL": 0, "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULTS": 0 }
        self.espreso.run(1, path, config)

    def espresoBinary(self, path, file):
        config = { "ENV::TESTING_LEVEL": 0, "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULTS": 0 }
        self.espreso.decompose(path)
        for input in glob.glob(path + "/ESDATA_*"):
            self.espreso.run(len(filter(os.path.isdir, glob.glob(input + "/*"))), path, config, [input])
            shutil.rmtree(input)

if __name__ == '__main__':

    openfoam, workbench, espresoBinary = TestCaseCreator.select(
        os.path.join(ROOT, "examples", "input", "openfoam"),
        os.path.join(ROOT, "examples", "input", "workbench"),
        os.path.join(ROOT, "examples", "input", "espresoBinaryFormat")
    )

    for subdirectory in openfoam:
        for name, path, file in TestCaseCreator.gather(subdirectory, ".ecf"):
            TestCaseCreator.create_test(ESPRESOInput, ESPRESOInput.openfoam, "openfoam_" + name, path, file)

    for subdirectory in workbench:
        for name, path, file in TestCaseCreator.gather(subdirectory, ".ecf"):
            TestCaseCreator.create_test(ESPRESOInput, ESPRESOInput.workbench, "workbench_" + name, path, file)

    for subdirectory in espresoBinary:
        for name, path, file in TestCaseCreator.gather(subdirectory, ".ecf", "decomposer.ecf"):
            TestCaseCreator.create_test(ESPRESOInput, ESPRESOInput.espresoBinary, "espresoBinary_" + name, path, file)

    unittest.main()
