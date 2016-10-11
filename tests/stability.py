
from utils import *
import unittest

class ESPRESOTests(unittest.TestCase):

    espreso = Espreso()
    cube = "tests/examples/linearElasticity/cube"

    def ordinary_check(self, config, info, oscilation=3):
        info.precision(1e-4)
        info.iterations(4)
        #info.oscilation(True, oscilation);

    def regular_cube(self, procs, config, args):
        config["ITERATIONS"] = 300
        config["EPSILON"] = 1e-4
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "regular_fixed_bottom.txt"
        info = RunInfo(self.espreso.output(procs, self.cube, config, args + [2, 2, 2, 2, 2, 2]))
        self.ordinary_check(config, info)

    def metis_cube(self, procs, config, args):
        config["ITERATIONS"] = 500
        config["EPSILON"] = 1e-4
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "metis_fixed_bottom.txt"
        info = RunInfo(self.espreso.output(procs, self.cube, config, args + [2, 2, 2, 5, 5, 5]))
        self.ordinary_check(config, info)

    def metis_cube_with_cyclic_edge(self, procs, config, args):
        config["ITERATIONS"] = 300
        config["EPSILON"] = 1e-4
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "metis_fixed_bottom.txt"
        info = RunInfo(self.espreso.output(procs, self.cube, config, args + [2, 1, 1, 3, 6, 6]))
        self.ordinary_check(config, info)

    def regular_esdata(self, procs, config, args):
        config["PATH"] = "regular_fixed_bottom.txt"
        dArgs = args[0:1] + [1, 1, 1] + args[1:4] + [int(12 / args[1]), int(12 / args[2]), int(12 / args[3]), "REGULAR_" ]
        self.espreso.decompose(self.cube, config, dArgs)

        config["ITERATIONS"] = 300
        config["EPSILON"] = 1e-4
        config["SUBDOMAINS"] = 5
        config["INPUT"] = "ESDATA"
        config["PATH"] = "REGULAR_" + str(procs)
        info = RunInfo(self.espreso.output(procs, self.cube, config, []))
        self.ordinary_check(config, info)

    def metis_esdata(self, procs, config, args):
        config["PATH"] = "metis_fixed_bottom.txt"
        dArgs = args[0:1] + [1, 1, 1] + args[1:4] + [int(12 / args[1]), int(12 / args[2]), int(12 / args[3]), "METIS_" ]
        self.espreso.decompose(self.cube, config, dArgs)

        config["ITERATIONS"] = 500
        config["EPSILON"] = 1e-4
        config["SUBDOMAINS"] = 4
        config["INPUT"] = "ESDATA"
        config["PATH"] = "METIS_" + str(procs)
        info = RunInfo(self.espreso.output(procs, self.cube, config, []))
        self.ordinary_check(config, info)

if __name__ == '__main__':

    def create_instance(config, example):
        if config["FETI_METHOD"] == "TOTAL_FETI" and config["B0_TYPE"] == "KERNELS":
            return
        procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
        args = [example["ETYPE"]] + example["CLUSTERS"]
        name = "_".join(str(x) for x in args + config.values())
        config["VERBOSE_LEVEL"] = 1
        config["TESTING_LEVEL"] = 1
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.regular_cube, name, procs, config, args)
        #TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.metis_cube, name + "_METIS", procs, config, args)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.metis_cube_with_cyclic_edge, name + "_METIS_TWO_SUBDOMAINS", procs, config, args)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.regular_esdata, name + "_ESDATA", procs, config, args)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.metis_esdata, name + "_METIS_ESDATA", procs, config, args)

    config = {
      "FETI_METHOD": [ "TOTAL_FETI", "HYBRID_FETI" ],
      "PRECONDITIONER": [ "NONE", "LUMPED", "WEIGHT_FUNCTION", "DIRICHLET" ],
      "REGULARIZATION": [ "FIX_POINTS", "NULL_PIVOTS" ],
      "B0_TYPE": [ "CORNERS", "KERNELS" ],
      "CGSOLVER": [ "STANDARD", "PIPELINED", "FULL_ORTOGONAL" ]
    }

    example = {
      "ETYPE":  [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5", "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ],
      "CLUSTERS": [ [1, 1, 1], [1, 1, 4], [2, 2, 2] ]
    }

    TestCaseCreator.iterate(create_instance, config, example)
    unittest.main()
