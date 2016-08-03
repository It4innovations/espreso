
from optparse import OptionParser
from utils import *
import unittest


class ESPRESOExhaustiveTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.espreso = Espreso("__ESPRESO__")
        cls.espreso.install(cls.config)

    @classmethod
    def tearDownClass(cls):
        cls.espreso.remove()

    def simple(self, esconfig):
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "metis_cube_elasticity_fixed_bottom.txt"
        self.espreso.run(4, "examples/meshgenerator/", esconfig, [0,  2, 2, 1,  2, 2, 2,  4, 4, 4])


if __name__ == '__main__':

    espreso_config = {
      "FETI_METHOD": [ "TOTAL_FETI", "HYBRID_FETI" ],
      "PRECONDITIONER": [ "LUMPED", "DIRICHLET" ]
    }

    build_config = {
      "SOLVER": [ "MKL", "PARDISO", "CUDA" ],
      "LIBTYPE": [ "SHARED", "STATIC" ],
      "INT_WIDTH": [ "32", "64" ]
    }

    testClasses = unittest.TestSuite()
    loader = unittest.TestLoader()
    opt = OptionParser()
    opt.add_option("-v", "--verbose", action="store_const", const=2, dest="verbose")


    def create_test(config):
        name = "_".join(config.values())
        TestCaseCreator.create_test(ESPRESOExhaustiveTests, ESPRESOExhaustiveTests.simple, name, config)

    def create_class(config):
        name = "_".join(config.values())
        testClass = type(name, (ESPRESOExhaustiveTests, ), {'config': config})
        testClasses.addTest(loader.loadTestsFromTestCase(testClass))


    """
    Iterate over ESPRESO configuration and build configuration
    """
    TestCaseCreator.iterate(create_test, espreso_config)
    TestCaseCreator.iterate(create_class, build_config)


    runner = unittest.TextTestRunner( verbosity = opt.parse_args()[0].verbose )
    runner.run( testClasses )