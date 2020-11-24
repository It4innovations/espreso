
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 3, 2, 1, 1, 2, 3, 4, 3, 4, "solver", "method" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA8", "HEXA20", "TETRA4", "TETRA10", "PRISMA6", "PRISMA15", "PYRAMID5", "PYRAMID13" ]:
        for method in [ "TOTAL_FETI", "HYBRID_FETI" ]:
            yield run, etype, "FETI", method
        if ESPRESOTest.fast: return

def run(etype, solver, method):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[10] = solver
    ESPRESOTest.args[11] = method
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))
