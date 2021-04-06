
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "TOTAL_FETI", "cgsolver", "preconditioner", "regularization" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for cgsolver in [ "PCG", "pipePCG", "orthogonalPCG", "GMRES", "BICGSTAB" ]:
        for preconditioner in [ "NONE", "LUMPED", "WEIGHT_FUNCTION", "DIRICHLET" ]:
            for regularization in [ "ANALYTIC", "ALGEBRAIC" ]:
                yield run, cgsolver, preconditioner, regularization

def run(cgsolver, preconditioner, regularization):
    ESPRESOTest.args[1] = cgsolver
    ESPRESOTest.args[2] = preconditioner
    ESPRESOTest.args[3] = regularization
    ESPRESOTest.run()
