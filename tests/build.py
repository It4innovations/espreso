
import shutil

from utils import *
import unittest

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ESPRESO_ROOT = os.path.dirname(ESPRESO_TESTS)
ESPRESO_LIBS = os.path.join(ESPRESO_ROOT, "libs/")
EXAMPLES = os.path.join(ESPRESO_TESTS, "examples")
TEST_DIR = "__espreso__"
TEST_LIBS = os.path.join(TEST_DIR, "libs/")

class TestBuild(unittest.TestCase):

    def setUp(self):
        shutil.rmtree(os.path.join(EXAMPLES, TEST_DIR), ignore_errors=True)
        subprocess.call(["git", "clone", "git@code.it4i.cz:mec059/espreso.git", TEST_DIR, "-q"], cwd=EXAMPLES)
        if os.path.isfile(ESPRESO_LIBS + "libpardiso500-INTEL120-X86-64.so"):
            subprocess.call(["mkdir", TEST_LIBS], cwd=EXAMPLES)
            subprocess.call(["cp", ESPRESO_LIBS + "libpardiso500-INTEL120-X86-64.so", TEST_LIBS], cwd=EXAMPLES)
            subprocess.call(["cp", ESPRESO_LIBS + "libifcore.a", TEST_LIBS], cwd=EXAMPLES)

    def tearDown(self):
        shutil.rmtree(os.path.join(EXAMPLES, TEST_DIR))

    def build(self, config):
        def check(result, error, method):
            success = False
            for line in result.splitlines():
                if line.find("'" + method + "' finished successfully") != -1:
                    success = True

            self.assertEqual(error, "", error)
            self.assertTrue(success, result)


        espreso = Espreso(TEST_DIR)
        result, error = espreso.configure(config)
        check(result, error, "configure")
        result, error = espreso.build()
        check(result, error, "install")

        esconfig = {
          "FETI_METHOD": "HYBRID_FETI",
          "PRECONDITIONER": "DIRICHLET"
        }

        if config["SOLVER"] == "CUDA":
            espreso.append_env("LD_LIBRARY_PATH", "/usr/local/cuda-7.5/lib64")

        result, error = espreso.run(4, "examples/meshgenerator/metis_cube_elasticity_fixed_bottom.txt", "GENERATOR", [ 0, 2, 2, 1, 2, 2, 2, 4, 4, 4], esconfig)
        self.assertEqual(result, "", result)
        self.assertEqual(error, "", error)


if __name__ == '__main__':

    def create_instance(settings):
        name = "_".join(settings.values())
        TestCaseCreator.create_test(TestBuild, TestBuild.build, name, settings.get())

    settings = Iterator({
      "SOLVER": [ "MKL", "CUDA" ],
      "LIBTYPE": [ "SHARED", "STATIC" ],
      "INT_WIDTH": [ "32", "64" ]
    })

    if os.path.isfile(ESPRESO_LIBS + "libpardiso500-INTEL120-X86-64.so"):
        settings.items["SOLVER"].append("PARDISO")


    TestCaseCreator.iterate(create_instance, settings)
    unittest.main()