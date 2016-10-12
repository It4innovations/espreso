
from utils import *
import unittest

class ESPRESOTests(unittest.TestCase):

    espreso = Espreso()
    cube = "tests/examples/linearElasticity/cube"

    def ordinary_check(self, config, info):
        if "ITERATIONS" in config and config["ITERATIONS"] > 3:
            info.iterations(3)
        if "EPSILON" in config:
            info.precision(config["EPSILON"])

    def regular_cube(self, procs, config, args):
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "regular_fixed_bottom.txt"
        info = RunInfo(self.espreso.output(procs, self.cube, config, args))
        self.ordinary_check(config, info)

    def metis_cube(self, procs, config, args):
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "metis_fixed_bottom.txt"
        info = RunInfo(self.espreso.output(procs, self.cube, config, args))
        self.ordinary_check(config, info)

    def esdata(self, procs, config, args):
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "metis_fixed_bottom.txt"
        self.espreso.decompose(self.cube, config, args + [ "METIS_" ])

        procs = reduce(lambda x, y: x * y, args[4:7])
        config["INPUT"] = "ESDATA"
        config["PATH"] = "METIS_" + str(procs)
        info = RunInfo(self.espreso.output(procs, self.cube, config, []))
        self.ordinary_check(config, info)

if __name__ == '__main__':

    def parameters(config, example):
        procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
        args = [example["ETYPE"]] + example["CLUSTERS"] + example["ARGS"]
        name = "_".join(str(x) for x in args + config.values())
        config["VERBOSE_LEVEL"] = 1
        config["TESTING_LEVEL"] = 1
        return name, procs, args

    def regular_cube(config, example):
        name, procs, args = parameters(config, example)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.regular_cube, "REGULAR_CUBE_" + name, procs, config, args)

    def metis_cube(config, example):
        name, procs, args = parameters(config, example)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.metis_cube, "METIS_CUBE_" + name, procs, config, args)

    def esdata(config, example):
        name, procs, args = parameters(config, example)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.esdata, "ESDATA_" + name, procs, config, args)

    FETI_METHODS = [ "TOTAL_FETI", "HYBRID_FETI" ]
    PRECONDITIONERS = [ "NONE", "LUMPED", "WEIGHT_FUNCTION", "DIRICHLET" ]
    REGULARIZATIONS = [ "FIX_POINTS", "NULL_PIVOTS" ]
    B0_TYPES = [ "CORNERS", "KERNELS" ]
    CGSOLVERS = [ "STANDARD", "PIPELINED", "FULL_ORTOGONAL" ]
    ETYPES = [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5", "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ]

    # Test generator for corners and fix points
    TestCaseCreator.iterate(
        regular_cube,
        {
            "FETI_METHOD": [ "HYBRID_FETI" ],
            "PRECONDITIONER": [ "DIRICHLET" ],
            "REGULARIZATION": [ "FIX_POINTS" ],
            "B0_TYPE": [ "CORNERS" ],
            "ITERATIONS": [ 40 ],
            "EPSILON": [ 1e-5 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 2, 4, 3, 5, 3, 3] ]
        }
    )

    # Dirichlet preconditioner for 1 subdomain has to compute results in 1 step
    TestCaseCreator.iterate(
        regular_cube,
        {
            "FETI_METHOD": [ "TOTAL_FETI" ],
            "PRECONDITIONER": [ "DIRICHLET" ],
            "REGULARIZATION": [ "FIX_POINTS" ],
            "ITERATIONS": [ 2 ],
            "EPSILON": [ 1e-5 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 1, 1, 1, 3, 4, 2] ]
        }
    )

    # Check computing corners on cyclic edges
    TestCaseCreator.iterate(
        metis_cube,
        {
            "FETI_METHOD": [ "HYBRID_FETI" ],
            "PRECONDITIONER": [ "DIRICHLET" ],
            "B0_TYPE": [ "CORNERS" ],
            "ITERATIONS": [ 30 ],
            "EPSILON": [ 1e-5 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 2, 1, 1, 2, 4, 4] ]
        }
    )

    # Check save/load ESDATA
    TestCaseCreator.iterate(
        esdata,
        {
            "FETI_METHOD": [ "TOTAL_FETI" ],
            "PRECONDITIONER": [ "DIRICHLET" ],
            "ITERATIONS": [ 30 ],
            "EPSILON": [ 1e-5 ],
            "SUBDOMAINS": [ 4 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 2, 2, 2, 6, 6, 6] ]
        }
    )

    # Main test for TOTAL FETI - loop over all settings
    TestCaseCreator.iterate(
        metis_cube,
        {
            "FETI_METHOD": [ "TOTAL_FETI" ],
            "PRECONDITIONER": PRECONDITIONERS,
            "REGULARIZATION": REGULARIZATIONS,
            "CGSOLVER": CGSOLVERS,
            "ITERATIONS": [ 600 ],
            "EPSILON": [ 1e-4 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [1, 2, 2] ],
            "ARGS": [ [3, 1, 2, 4, 6, 3] ]
        }
    )

    # Main test for HYBRID FETI - loop over all settings
    TestCaseCreator.iterate(
        metis_cube,
        {
            "FETI_METHOD": [ "HYBRID_FETI" ],
            "PRECONDITIONER": PRECONDITIONERS,
            "REGULARIZATION": REGULARIZATIONS,
            "B0_TYPE": B0_TYPES,
            "CGSOLVER": CGSOLVERS,
            "ITERATIONS": [ 600 ],
            "EPSILON": [ 1e-4 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [1, 2, 2] ],
            "ARGS": [ [3, 1, 2, 4, 6, 3] ]
        }
    )

    unittest.main()
