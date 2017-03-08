
from utils import *
import unittest

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(ESPRESO_TESTS)

class ESPRESOBenchmarks(unittest.TestCase):

    espreso = Espreso()

    def benchmark(self, path, file):
        config = { "ENV::TESTING_LEVEL": 1, "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULTS": 0 }
        for line in [ line.rstrip('\n') for line in open(os.path.join(path, file)) ]:
            param, value = line.split("=")
            if param.strip() == "PROCS":
                procs = int(value)
            elif param.strip() == "ARGS":
                args = value.split()
            else:
                config["RESULTS::" + param.strip()] = value.strip()

        self.espreso.run(procs, path, config, args)


if __name__ == '__main__':

    benchmarks = TestCaseCreator.select(os.path.join(ROOT, "benchmarks"))

    for subdirectory in benchmarks:
        for name, path, file in TestCaseCreator.gather(subdirectory, ".test"):
            TestCaseCreator.create_test(ESPRESOBenchmarks, ESPRESOBenchmarks.benchmark, name, path, file)

    unittest.main()
