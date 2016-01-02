
from utils import *
import unittest

class TestMeshMethods(unittest.TestCase):

    def test_simple(self):
        result = Espreso(CubeGenerator(0, (1, 1, 1), (3, 3, 3)),
                Assembler("FEM", "LinearElasticity"),
                Solver(100, 0.001)
                ).run(2)

        # check result

if __name__ == '__main__':
    unittest.main()