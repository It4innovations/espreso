
from utils import *
import unittest

class TestESPRESO(unittest.TestCase):

    def test_simple(self):
        result, error = Espreso("linearElasticity/cube").run(1, [0, 1, 1, 1, 2, 2, 2, 3, 3, 3])
        self.assertEqual(error, "", "ESPRESO failed computation with the following error: \n" + error)

if __name__ == '__main__':
    unittest.main()