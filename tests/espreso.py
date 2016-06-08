
from utils import *
import unittest

class TestESPRESO(unittest.TestCase):

    def stability(self, procs, args, config):
        result, error = Espreso("linearElasticity/cube").run(procs, "regular_fixed_bottom.txt", "GENERATOR", args, config)
        self.assertEqual(result, "", result)
        self.assertEqual(error, "", error)


if __name__ == '__main__':
    def create_instance(example, config):
            procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
            args = [example["ETYPE"]] + example["CLUSTERS"]
            name = "_".join(str(x) for x in args + config.values())
            TestCaseCreator.create_test(TestESPRESO, TestESPRESO.stability, name, procs, args, config.get())

    example = Iterator({
      "ETYPE":  [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5" ],
      "CLUSTERS": [ [1, 1, 1], [1, 1, 4], [2, 2, 2] ]
    })

    config = Iterator({
      "FETI_METHOD": [ "0", "1" ],
      "PRECONDITIONER": [ "0", "1", "2", "3" ],
      "REGULARIZATION": [ "0", "1" ]
    })

    TestCaseCreator.iterate(create_instance, example, config)
    unittest.main()