
from utils import *
import unittest

class TestESPRESO(unittest.TestCase):

    def test_simple(self):
        result = Espreso("linearElasticity/cube").run(1, [0, 1, 1, 1, 2, 2, 2, 3, 3, 3, "-v"])

        # check result

if __name__ == '__main__':
    unittest.main()