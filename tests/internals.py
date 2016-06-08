
from utils import *
import unittest

class TestESPRESO(unittest.TestCase):

    def test_parameters(self):
        result, error = Espreso("parameters").run(1, "correct.txt", "GENERATOR", [], {})
        self.assertEqual(result, "", result)
        self.assertEqual(error, "", error)

        for file in os.listdir(os.path.join(EXAMPLES, "parameters")):
            if file.startswith("error"):
                result, error = Espreso("parameters").run(1, file, "GENERATOR", [], {})
                self.assertEqual(result, "", result)
                self.assertNotEqual(error, "", error)

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
            result, error = Espreso("parameters").run(1, "correct.txt", "GENERATOR", [], { key: value })
            self.assertEqual(result, "", result)
            self.assertEqual(error, "", error)

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
            result, error = Espreso("parameters").run(1, "correct.txt", "GENERATOR", [], {key: value})
            self.assertEqual(result, "", result)
            self.assertEqual(error, "Parameter '" + key + "' has a wrong value '" + value + "'.\n")

        parameters = {
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

        result, error = Espreso("parameters", ).run(1, "correct.txt", "GENERATOR", [], parameters)
        self.assertEqual(error, "", error)
        found = 0
        for line in result.splitlines():
            for out in parameters_output:
                if line.find(out) != -1:
                    found += 1

        self.assertEqual(found, len(parameters_output), "Not all parameters are parsed correctly.")



if __name__ == '__main__':
    unittest.main()