
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
            TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.stability, name, procs, config.get(), args)

    config = Iterator({
      "FETI_METHOD": [ "0", "1" ],
      "PRECONDITIONER": [ "0", "1", "2", "3" ],
      "REGULARIZATION": [ "0", "1" ],
    })

    example = Iterator({
      "ETYPE":  [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5" ],
      "CLUSTERS": [ [1, 1, 1], [1, 1, 4], [2, 2, 2] ]
    })

    TestCaseCreator.iterate(create_instance, config, example)
    unittest.main()