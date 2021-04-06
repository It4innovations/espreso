
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "method", "cgsolver", "preconditioner", "ALGEBRAIC", "B0 type" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for cgsolver in [ "GMRES", "BICGSTAB" ]:
        for preconditioner in [ "NONE", "LUMPED", "WEIGHT_FUNCTION", "DIRICHLET" ]:
            yield run, "TOTAL_FETI", cgsolver, preconditioner, "KERNELS"
            for B0_type in [ "CORNERS", "KERNELS" ]:
                yield run, "HYBRID_FETI", cgsolver, preconditioner, B0_type

def run(method, cgsolver, preconditioner, B0_type):
    ESPRESOTest.args[0] = method
    ESPRESOTest.args[1] = cgsolver
    ESPRESOTest.args[2] = preconditioner
    ESPRESOTest.args[4] = B0_type
    ESPRESOTest.run()
