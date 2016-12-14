
from utils import *
import unittest

class ESPRESOTests(unittest.TestCase):

    espreso = Espreso()
    regular_cube  = "tests/examples/regular_cube"
    iregular_cube = "tests/examples/iregular_cube"

    def ordinary_check(self, config, info):
        if "ITERATIONS" in config and config["ITERATIONS"] > 3:
            info.iterations(3)
        if "EPSILON" in config:
            info.precision(config["EPSILON"])

    def cube(self, example, procs, config, args):
        info = RunInfo(self.espreso.output(procs, example, config, args))
        self.ordinary_check(config, info)

    def esdata(self, procs, config, args):
        self.espreso.decompose(self.regular_cube, config, args + [ "METIS_" ])

        print "decomposed"
        procs = reduce(lambda x, y: x * y, args[4:7])
        config["INPUT"] = "ESDATA"
        config["ESDATA::PATH"] = "DECOMPOSITION" + str(procs)
        info = RunInfo(self.espreso.output(procs, self.regular_cube, config, args))
        self.ordinary_check(config, info)

if __name__ == '__main__':

    def parameters(config, example):
        procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
        args = [example["ETYPE"]] + example["CLUSTERS"] + example["ARGS"]
        name = "_".join(str(x) for x in args + config.values())
        config["ENV::VERBOSE_LEVEL"] = 1
        config["ENV::TESTING_LEVEL"] = 1
        return name, procs, args

    def regular_cube(config, example):
        name, procs, args = parameters(config, example)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.cube, "REGULAR_CUBE_" + name, ESPRESOTests.regular_cube, procs, config, args)

    def metis_cube(config, example):
        name, procs, args = parameters(config, example)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.cube, "METIS_CUBE_" + name, ESPRESOTests.iregular_cube, procs, config, args)

    def esdata(config, example):
        name, procs, args = parameters(config, example)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.esdata, "ESDATA_" + name, procs, config, args)

    FETI_METHODS = [ "TOTAL_FETI", "HYBRID_FETI" ]
    PRECONDITIONERS = [ "NONE", "LUMPED", "WEIGHT_FUNCTION", "DIRICHLET" ]
    REGULARIZATIONS = [ "FIX_POINTS", "NULL_PIVOTS" ]
    B0_TYPES = [ "CORNERS", "KERNELS" ]
    CGSOLVERS = [ "PCG", "pipePCG", "orthogonalPCG" ]
    ETYPES = [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5", "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ]

    SOLVER = "LINEAR_ELASTICITY_3D::ESPRESO::"

    # Test generator for corners and fix points
    TestCaseCreator.iterate(
        regular_cube,
        {
            SOLVER + "METHOD": [ "HYBRID_FETI" ],
            SOLVER + "PRECONDITIONER": [ "DIRICHLET" ],
            SOLVER + "REGULARIZATION": [ "FIX_POINTS" ],
            SOLVER + "B0_TYPE": [ "CORNERS" ],
            SOLVER + "ITERATIONS": [ 40 ],
            SOLVER + "EPSILON": [ 1e-5 ]
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
            SOLVER + "METHOD": [ "TOTAL_FETI" ],
            SOLVER + "PRECONDITIONER": [ "DIRICHLET" ],
            SOLVER + "REGULARIZATION": [ "FIX_POINTS" ],
            SOLVER + "ITERATIONS": [ 2 ],
            SOLVER + "EPSILON": [ 1e-5 ]
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
            SOLVER + "METHOD": [ "HYBRID_FETI" ],
            SOLVER + "PRECONDITIONER": [ "DIRICHLET" ],
            SOLVER + "B0_TYPE": [ "CORNERS" ],
            SOLVER + "ITERATIONS": [ 30 ],
            SOLVER + "EPSILON": [ 1e-5 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 2, 1, 1, 2, 4, 4] ]
        }
    )

# TODO: FIX ESDATA
#    # Check save/load ESDATA
#    TestCaseCreator.iterate(
#        esdata,
#        {
#            SOLVER + "METHOD": [ "TOTAL_FETI" ],
#            SOLVER + "PRECONDITIONER": [ "DIRICHLET" ],
#            SOLVER + "ITERATIONS": [ 30 ],
#            SOLVER + "EPSILON": [ 1e-5 ],
#            "ESDATA::DOMAINS": [ 4 ]
#        },
#        {
#            "ETYPE": [ "HEXA8" ], #ETYPES,
#            "CLUSTERS": [ [1, 1, 1] ],
#            "ARGS": [ [ 2, 2, 2, 6, 6, 6] ]
#        }
#    )

    # Main test for TOTAL FETI - loop over all settings
    TestCaseCreator.iterate(
        metis_cube,
        {
            SOLVER + "METHOD": [ "TOTAL_FETI" ],
            SOLVER + "PRECONDITIONER": PRECONDITIONERS,
            SOLVER + "REGULARIZATION": REGULARIZATIONS,
            SOLVER + "SOLVER": CGSOLVERS,
            SOLVER + "ITERATIONS": [ 600 ],
            SOLVER + "EPSILON": [ 1e-4 ]
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
            SOLVER + "METHOD": [ "HYBRID_FETI" ],
            SOLVER + "PRECONDITIONER": PRECONDITIONERS,
            SOLVER + "REGULARIZATION": REGULARIZATIONS,
            SOLVER + "B0_TYPE": B0_TYPES,
            SOLVER + "SOLVER": CGSOLVERS,
            SOLVER + "ITERATIONS": [ 600 ],
            SOLVER + "EPSILON": [ 1e-4 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [1, 2, 2] ],
            "ARGS": [ [3, 1, 2, 4, 6, 3] ]
        }
    )

    unittest.main()
