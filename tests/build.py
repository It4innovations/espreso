
import shutil

from utils import *
import unittest

class ESPRESOBuildTests(unittest.TestCase):

    espreso = Espreso("__ESPRESO__")

    @classmethod
    def tearDownClass(cls):
        cls.espreso.remove()

    def setUp(self):
        self.espreso.waf(["distclean"])

    def build(self, config):
        self.espreso.install(config)

        esconfig = {
          "FETI_METHOD": "HYBRID_FETI",
          "PRECONDITIONER": "DIRICHLET",
          "INPUT": "GENERATOR",
          "PATH": "metis_cube_elasticity_fixed_bottom.txt"
        }

        self.espreso.run(4, "examples/meshgenerator/", esconfig, [0,  2, 2, 1,  2, 2, 2,  4, 4, 4])


if __name__ == '__main__':

    def create_instance(settings):
        name = "_".join(settings.values())
        TestCaseCreator.create_test(ESPRESOBuildTests, ESPRESOBuildTests.build, name, settings)

    settings = {
      "SOLVER": [ "MKL", "CUDA", "PARDISO" ],
      "LIBTYPE": [ "SHARED", "STATIC" ],
      "INT_WIDTH": [ "32", "64" ]
    }


    TestCaseCreator.iterate(create_instance, settings)
    unittest.main()