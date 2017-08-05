
import os
import sys
import shutil

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(ESPRESO_TESTS)
sys.path.append(os.path.join(ESPRESO_TESTS, "utils"))

from testing import *
import unittest

class ESPRESOBenchmarks(unittest.TestCase):

    espreso = Espreso()

    def benchmark(self, path, file):
        self.espreso.run(
            self.espreso.get_processes(os.path.join(path, file)),
            path,
            { "config": file, "ENV::TESTING_LEVEL": 0, "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULTS": 0 }
        )

        self.espreso.compare_monitors(
            os.path.join(path, file.replace(".ecf", ".emr")),
            os.path.join(path, "results", "last", file.replace(".ecf", ".emr"))
        )


if __name__ == '__main__':

    benchmarks = TestCaseCreator.select(os.path.join(ROOT, "benchmarks"))

    for subdirectory in benchmarks:
        for name, path, file in TestCaseCreator.gather(subdirectory, ".ecf"):
            TestCaseCreator.create_test(ESPRESOBenchmarks, ESPRESOBenchmarks.benchmark, name, path, file)

    unittest.main()
