
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "TOTAL_FETI", "cgsolver", "preconditioner", "ALGEBRAIC" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for cgsolver in [ "GMRES", "BICGSTAB" ]:
        for preconditioner in [ "NONE", "LUMPED", "WEIGHT_FUNCTION", "DIRICHLET" ]:
            yield run, cgsolver, preconditioner

def run(cgsolver, preconditioner):
    ESPRESOTest.args[1] = cgsolver
    ESPRESOTest.args[2] = preconditioner
    ESPRESOTest.run()
