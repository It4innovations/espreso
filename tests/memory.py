
from utils import *
import unittest

class ESPRESOTests(unittest.TestCase):

    espreso = Espreso()
    cube = "tests/examples/linearElasticity/cube"

    def valgrind(self, procs, config, args):
        config["ITERATIONS"] = 50
        config["EPSILON"] = 1e-2
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "metis_fixed_bottom.txt"
        config["VERBOSE_LEVEL"] = 0
        config["TESTING_LEVEL"] = 0
        self.espreso.valgrind(procs, self.cube, config, args + [2, 2, 1, 3, 3, 5])

if __name__ == '__main__':

    def create_instance(config, example):
        if config["FETI_METHOD"] == "TOTAL_FETI" and config["B0_TYPE"] == "KERNELS":
            return
        procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
        args = [example["ETYPE"]] + example["CLUSTERS"]
        name = "_".join(str(x) for x in args + config.values())
        config["VERBOSE_LEVEL"] = 1
        config["TESTING_LEVEL"] = 1
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.valgrind, name + "_VALGRIND", procs, config, args)

    config = {
      "FETI_METHOD": [ "TOTAL_FETI", "HYBRID_FETI" ],
      "PRECONDITIONER": [ "NONE", "LUMPED", "WEIGHT_FUNCTION", "DIRICHLET" ],
      "REGULARIZATION": [ "FIX_POINTS", "NULL_PIVOTS" ],
      "B0_TYPE": [ "CORNERS", "KERNELS" ],
      "CGSOLVER": [ "STANDARD", "PIPELINED", "FULL_ORTOGONAL" ]
    }

    TestCaseCreator.iterate(create_instance, config, { "ETYPE": [ "HEXA8" ], "CLUSTERS": [ [1, 2, 2] ] })

    config = {
      "FETI_METHOD": [ "TOTAL_FETI", "HYBRID_FETI" ],
      "PRECONDITIONER": [ "DIRICHLET" ],
      "REGULARIZATION": [ "FIX_POINTS", "NULL_PIVOTS" ],
      "B0_TYPE": [ "CORNERS", "KERNELS" ],
      "CGSOLVER": [ "STANDARD" ]
    }

    example = {
      "ETYPE":  [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5", "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ],
      "CLUSTERS": [ [1, 2, 2] ]
    }
    TestCaseCreator.iterate(create_instance, config, example)

    unittest.main()
