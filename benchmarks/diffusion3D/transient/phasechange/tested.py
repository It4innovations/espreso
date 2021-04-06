
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 4, 1, 1, 3, 1, 1, 4, 4, 4, "solver" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ]:
        for solver in [ "FETI", "HYPRE", "MKLPDSS" ]:
            yield run, etype, solver
        if ESPRESOTest.fast: return

def run(etype, solver):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[10] = solver
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
