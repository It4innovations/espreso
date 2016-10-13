
from utils import *
import unittest
import os
import glob

class ESPRESOTests(unittest.TestCase):

    espreso = Espreso()
    cube = "tests/examples/linearElasticity/cube"

    def regular_cube(self, procs, config, args):
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "regular_fixed_bottom.txt"
        self.espreso.valgrind(procs, self.cube, config, args)


if __name__ == '__main__':

    def parameters(config, example):
        procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
        args = [example["ETYPE"]] + example["CLUSTERS"] + example["ARGS"]
        name = "_".join(str(x) for x in args + config.values())
        config["SAVE_RESULTS"] = 1
        return name, procs, args

    def regular_cube(config, example):
        name, procs, args = parameters(config, example)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.regular_cube, "VALGRIND_OUTPUT" + name, procs, config, args)

    OUTPUT_FORMATS = [ "VTK_LEGACY_FORMAT", "VTK_BINARY_FORMAT", "VTK_MULTIBLOCK_FORMAT", "ENSIGHT_FORMAT" ]

    # Test output format
    TestCaseCreator.iterate(
        regular_cube,
        {
            "OUTPUT_FORMAT": OUTPUT_FORMATS,
            "OUTPUT_COMPRESSION": [ 0, 1 ],
            "OUTPUT_DECIMATION": [ 0, 0.5 ]
        },
        {
            "ETYPE": [ "HEXA8" ],
            "CLUSTERS": [ [2, 1, 1] ],
            "ARGS": [ [ 2, 1, 1, 1, 2, 2] ]
        }
    )

    unittest.main()
