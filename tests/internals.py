
from utils import *
import unittest

class ESPRESOInternalTests(unittest.TestCase):

    espreso = Espreso()

    def test_parameters(self):
        esconfig = {
            "INPUT": "GENERATOR",
            "PATH": "correct.txt",
        }
        example_dir = "tests/examples/parameters"

        self.espreso.run(1, example_dir, esconfig)

        for file in os.listdir(os.path.dirname(os.path.abspath(__file__)) + "/examples/parameters"):
            if file.startswith("error"):
                esconfig["PATH"] = file
                self.espreso.fail(1, example_dir, esconfig)

        esconfig["PATH"] = "correct.txt"

        correct = [
            ("INPUT", "4"),
            ("INPUT", "GENERATOR"),
            ("INPUT", "GENErator"),
            ("INPUT", "gENeRAtOR"),
            ("REDUNDANT_LAGRANGE", "0"),
            ("REDUNDANT_LAGRANGE", "1"),
            ("SUBDOMAINS", "10"),
            ("EPSILON", "0.005"),
            ("EPSILON", "1e-3"),
            ("PRECONDITIONER", "DIRICHLET"),
            ("PRECONDITIONER", "3"),
            ("PRECONDITIONER", "DIRICHlet"),
            ("PRECONDITIONER", "dIRICHLET"),
        ]

        for key, value in correct:
            config = esconfig.copy()
            config.update({key: value})
            self.espreso.run(1, example_dir, config)

        incorrect = [
            ("INPUT", "5"),
            ("INPUT", "GENERATO"),
            ("INPUT", "GENErators"),
            ("INPUT", "gENeRAtOR4"),
            ("INPUT", "4gENeRAtOR"),
            ("INPUT", "-1"),
            ("INPUT", "0.1"),
            ("INPUT", "4.0"),
            ("INPUT", "4ff"),
            ("INPUT", "f4f"),
            ("INPUT", "4-4"),
            ("SUBDOMAINS", "10.0"),
            ("EPSILON", "Hello"),
            ("EPSILON", "ee4"),
            ("REDUNDANT_LAGRANGE", "True"),
            ("REDUNDANT_LAGRANGE", "False"),
            ("USE_SCHUR_COMPLEMENT", "x"),
        ]

        for key, value in incorrect:
            config = esconfig.copy()
            config.update({key: value})
            self.espreso.fail(1, example_dir, config)

        parameters = {
            "INPUT": "GENERATOR",
            "PATH": "correct.txt",
            "SUBDOMAINS": "1024",
            "EPSILON": "1e-3",
            "USE_SCHUR_COMPLEMENT": "",
            "VERTEX_CORNERS": "0",
            "FACE_CORNERS": "1",
            "EDGE_CORNERS": "10",
            "VERBOSE_LEVEL": "3",
            "OUTPUT": "out",
        }
        parameters_output = [
            "SUBDOMAINS == 1024",
            "EPSILON == 0.001",
            "USE_SCHUR_COMPLEMENT == 1",
            "VERTEX_CORNERS == 0",
            "FACE_CORNERS == 1",
            "EDGE_CORNERS == 1",
            "VERBOSE_LEVEL == 3",
            "OUTPUT == out",
        ]

        output = self.espreso.output(1, example_dir, parameters)
        found = 0
        for line in output.splitlines():
            for out in parameters_output:
                if line.find(out) != -1:
                    found += 1

        self.assertEqual(found, len(parameters_output), "Not all parameters are parsed correctly.")



if __name__ == '__main__':
    unittest.main()