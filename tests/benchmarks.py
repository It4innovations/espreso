
from utils import *
import unittest
import glob
import string

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(ESPRESO_TESTS)

class ESPRESOBenchmarks(unittest.TestCase):

    espreso = Espreso()

    def benchmark(self, path):
        config = { "TESTING_LEVEL": 1, "VERBOSE_LEVEL": 0, "MEASURE_LEVEL": 0, "SAVE_RESULT": 0 }
        for test in  glob.glob(path + "/*.test"):
            for line in [ string.capwords(line.rstrip('\n')) for line in open(test) ]:
                param, value = line.split("=")
                if param.strip() == "Procs":
                    procs = int(value)
                if param.strip() == "Args":
                    args = [ int(i) for i in value.split() ]
                if param.strip() == "Norm":
                    config["NORM"] = float(value)

            self.espreso.run(procs, path, config, args)


if __name__ == '__main__':

    benchmarks = []
    for root, subFolders, files in os.walk(os.path.join(ROOT, "benchmarks")):
        if "espreso.test" in files:
            benchmarks.append(root)

    benchmarks.sort()
    for benchmark in benchmarks:
        TestCaseCreator.create_test(ESPRESOBenchmarks, ESPRESOBenchmarks.benchmark, os.path.basename(benchmark), benchmark)

    unittest.main()
