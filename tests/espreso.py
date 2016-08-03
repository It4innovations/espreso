
from utils import *
import unittest

class ESPRESOTests(unittest.TestCase):

    espreso = Espreso()

    def stability(self, procs, config, args):
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "regular_fixed_bottom.txt"
        self.espreso.run(procs, "tests/examples/linearElasticity/cube", config, args)


if __name__ == '__main__':

    def create_instance(config, example):
            procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
            args = [example["ETYPE"]] + example["CLUSTERS"]
            name = "_".join(str(x) for x in args + config.values())
            TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.stability, name, procs, config, args)

    config = {
      "FETI_METHOD": [ "TOTAL_FETI", "HYBRID_FETI" ],
      "PRECONDITIONER": [ "NONE", "LUMPED", "WEIGHT_FUNCTION", "DIRICHLET" ],
      "REGULARIZATION": [ "FIX_POINTS", "NULL_PIVOTS" ],
    }

    example = {
      "ETYPE":  [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5", "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ],
      "CLUSTERS": [ [1, 1, 1], [1, 1, 4], [2, 2, 2] ]
    }

    TestCaseCreator.iterate(create_instance, config, example)
    unittest.main()
