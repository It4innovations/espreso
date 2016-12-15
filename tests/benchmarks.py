
from utils import *
import unittest
import glob
import string

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(ESPRESO_TESTS)

class ESPRESOBenchmarks(unittest.TestCase):

    espreso = Espreso()

    def benchmark(self, path):
        config = { "ENV::TESTING_LEVEL": 3, "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULT": 0 }
        for test in  glob.glob(path + "/*.test"):
            for line in [ line.rstrip('\n') for line in open(test) ]:
                param, value = line.split("=")
                if param.strip() == "PROCS":
                    procs = int(value)
                elif param.strip() == "ARGS":
                    args = value.split()
                else:
                    config["RESULTS::" + param.strip()] = value.strip()

            self.espreso.run(procs, path, config, args)


if __name__ == '__main__':

    benchmarks = []
    for root, subFolders, files in os.walk(os.path.join(ROOT, "benchmarks")):
        for file in files:
            if file.endswith(".test"):
                benchmarks.append(root)
                break

    benchmarks.sort()
    for benchmark in benchmarks:
        TestCaseCreator.create_test(ESPRESOBenchmarks, ESPRESOBenchmarks.benchmark, os.path.basename(benchmark), benchmark)

    unittest.main()
