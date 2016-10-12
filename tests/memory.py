
from utils import *
import unittest

class ESPRESOTests(unittest.TestCase):

    espreso = Espreso()
    cube = "tests/examples/linearElasticity/cube"

    def valgrind(self, procs, config, args):
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "metis_fixed_bottom.txt"
        self.espreso.valgrind(procs, self.cube, config, args)

if __name__ == '__main__':

    def parameters(config, example):
        procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
        args = [example["ETYPE"]] + example["CLUSTERS"] + example["ARGS"]
        name = "_".join(str(x) for x in args + config.values())
        return name, procs, args

    def valgrind(config, example):
        name, procs, args = parameters(config, example)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.valgrind, "VALGRIND_" + name, procs, config, args)


    FETI_METHODS = [ "TOTAL_FETI", "HYBRID_FETI" ]
    PRECONDITIONERS = [ "NONE", "LUMPED", "WEIGHT_FUNCTION", "DIRICHLET" ]
    REGULARIZATIONS = [ "FIX_POINTS", "NULL_PIVOTS" ]
    B0_TYPES = [ "CORNERS", "KERNELS" ]
    CGSOLVERS = [ "STANDARD", "PIPELINED", "FULL_ORTOGONAL" ]
    ETYPES = [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5", "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ]

    # Test leaks of various preconditioners
    TestCaseCreator.iterate(
        valgrind,
        {
            "FETI_METHOD": [ "TOTAL_FETI" ],
            "PRECONDITIONER": PRECONDITIONERS,
            "ITERATIONS": [ 3 ]
        },
        {
            "ETYPE": [ "HEXA8" ],
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 2, 1, 1, 2, 4, 4] ]
        }
    )

    # Test leaks of various regularizations
    TestCaseCreator.iterate(
        valgrind,
        {
            "FETI_METHOD": [ "TOTAL_FETI" ],
            "REGULARIZATION": REGULARIZATIONS,
            "ITERATIONS": [ 3 ]
        },
        {
            "ETYPE": [ "HEXA8" ],
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 2, 1, 1, 2, 4, 4] ]
        }
    )

    # Test leaks of various cg solvers
    TestCaseCreator.iterate(
        valgrind,
        {
            "FETI_METHOD": [ "TOTAL_FETI" ],
            "CGSOLVER": CGSOLVERS,
            "ITERATIONS": [ 3 ]
        },
        {
            "ETYPE": [ "HEXA8" ],
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 2, 1, 1, 2, 4, 4] ]
        }
    )


    # Test leaks of various HYBRID FETI OBJECTS
    TestCaseCreator.iterate(
        valgrind,
        {
            "FETI_METHOD": [ "HYBRID_FETI" ],
            "PRECONDITIONER": [ "DIRICHLET" ],
            "B0_TYPE": B0_TYPES,
            "ITERATIONS": [ 3 ]
        },
        {
            "ETYPE": [ "HEXA8" ],
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 2, 1, 1, 2, 4, 4] ]
        }
    )

    # Test leaks of various elements
    TestCaseCreator.iterate(
        valgrind,
        {
            "FETI_METHOD": [ "HYBRID_FETI" ],
            "PRECONDITIONER": [ "DIRICHLET" ],
            "B0_TYPE": B0_TYPES,
            "ITERATIONS": [ 3 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [1, 1, 1] ],
            "ARGS": [ [ 2, 1, 1, 2, 4, 4] ]
        }
    )

    unittest.main()
