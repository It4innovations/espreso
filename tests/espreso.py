
from utils import *
import unittest

class TestESPRESO(unittest.TestCase):

    def stability(self, procs, args, config):
        def composeErrMsg(error):
            errorMsg = "ESPRESO fail computation for the following setting: \n"
            return errorMsg + error

        result, error = Espreso("linearElasticity/cube", "GENERATOR", "regular_fixed_bottom.txt", config).run(procs, args)
        self.assertEqual(result, "", result)
        self.assertEqual(error, "", error)

class TestCaseCreator:

    @staticmethod
    def createX(function, procs, args, config):
        def test_method(self):
            function(self, procs, args, config)

        name = "_".join(str(x) for x in args + config.values())
        setattr(TestESPRESO, 'test_' + name, test_method)
        test_method.__name__ = 'test_' + name

    @staticmethod
    def stability():
        example = Iterator({
          "ETYPE":  [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5" ],
          "CLUSTERS": [ [1, 1, 1], [1, 1, 4], [2, 2, 2] ]
        })

        config = Iterator({
          "FETI_METHOD": [ "0", "1" ],
          "PRECONDITIONER": [ "0", "1", "2", "3" ],
          "REGULARIZATION": [ "0", "1" ]
        })

        next_example = True
        while next_example:
            next_config = True
            while next_config:

                procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
                args = [example["ETYPE"]] + example["CLUSTERS"]
                TestCaseCreator.createX(TestESPRESO.stability, procs, args, config.get())

                next_config = config.next()
            next_example = example.next()


if __name__ == '__main__':
    TestCaseCreator.stability()
    unittest.main()