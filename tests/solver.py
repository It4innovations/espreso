
from utils import *
import unittest

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(ESPRESO_TESTS)

class ESPRESOSolver(unittest.TestCase):

    espreso = Espreso()

    def solver(self, path, file, args):
        config = { "ENV::TESTING_LEVEL": 1, "ENV::VERBOSE_LEVEL": 0, "ENV::MEASURE_LEVEL": 0, "OUTPUT::RESULTS": 0 }

        procs = 1
        for line in open(os.path.join(path, "espreso.ecf")):
            if line.strip().startswith("CLUSTER"):
                procs = procs * int(line.split()[1].strip(" ;"))

        self.espreso.run(procs, path, config, args)

if __name__ == '__main__':

    test_cases = os.path.join(ESPRESO_TESTS, "examples", "solver")

    for name, path, file in TestCaseCreator.gather(test_cases, ".test"):
        def create_test(*args):
            args = [ arg.values()[0] for arg in args ]
            suffix = os.path.relpath(path, test_cases).replace('/', '_')
            TestCaseCreator.create_test(ESPRESOSolver, ESPRESOSolver.solver, "_".join([suffix] + args), path, file, args)

        ranges = []
        for parameter, values in [ line.split("=") for line in open(os.path.join(path, file)) ]:
            ranges.append({ parameter.strip(): values.split() })

        TestCaseCreator.iterate(create_test, *ranges)

    unittest.main()
