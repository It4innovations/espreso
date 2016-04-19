
from utils import *
import unittest

class TestESPRESO(unittest.TestCase):

    def stability(self, id, procs, args, config):
        def composeErrMsg(error):
            errorMsg = "ESPRESO fail computation for the following setting: \n"
            return errorMsg + error

        result, error = Espreso("linearElasticity/cube", config).run(procs, args)
        self.assertEqual(result, "", result)
        self.assertEqual(error, "", error)


def _add_test(name, param1):
    def test_method(self):
        self.stability()

    setattr(TestESPRESO, 'test_' + name, test_method)
    test_method.__name__ = 'test_' + name

class TestCaseComposer:

    @staticmethod
    def stability():
        def create(name, id ,procs, args, config):
            def test_method(self):
                self.stability(id, procs, args, config)

            setattr(TestESPRESO, 'test_' + name, test_method)
            test_method.__name__ = 'test_' + name

        config = Configuration("stability");

        element_types = [ "HEXA8" ] #, "HEXA20", "TETRA4", "TETRA10", "PRISMA6", "PRISMA15", "PYRAMID5", "PYRAMID13" ]
        subdomains = [ 2, 2, 2 ]
        elements = [ 3, 3, 3 ]

        feti_method = [ "TFETI", "HTFETI" ]
        preconditioner = [ "NOPREC", "LUMPED", "WEIGHT", "DIRICHLET" ]
        regularization = [ "FIXPOINTS", "NULLPIVOTS" ]

        for method in range(0, len(feti_method)):
            for prec in range(0, len(preconditioner)):
                for reg in range(0, len(regularization)):

                    config.set("FETI_METHOD", method)
                    config.set("PRECONDITIONER", prec)
                    config.set("REGULARIZATION", reg)
                    suffix = "_".join([ feti_method[method], preconditioner[prec], regularization[reg] ]);
                    cfile = config.save("linearElasticity/cube", suffix)

                    for clusters in [ [1, 1, 1], [2, 1, 1], [1, 1, 4], [1, 8, 1], [2, 2, 2] ]:
                        for type in range(0, len(element_types)):
                            name = element_types[type] + "_args_" + "_".join(str(x) for x in clusters + subdomains + elements) + "_conf_" + suffix
                            procs = reduce(lambda x, y: x * y, clusters)
                            create(name, 1, procs, [ type ] + clusters + subdomains + elements, cfile)


if __name__ == '__main__':
    TestCaseComposer.stability()
    unittest.main()