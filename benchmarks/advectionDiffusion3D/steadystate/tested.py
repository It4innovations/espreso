
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 1, 4, 1, 2, 2, 2, 4, 6, 4, "solver" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5" ]:
        for solver in [ "FETI", "HYPRE", "MKLPDSS" ]:
            yield run, etype, solver
        if ESPRESOTest.fast: return

def run(etype, solver):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[10] = solver
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
