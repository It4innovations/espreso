
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.processes = 2
    ESPRESOTest.args = [ "etype", 2, 1, 8, 2, 10, 10, "solver" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "SQUARE8", "TRIANGLE6" ]:
        for solver in [ "FETI", "HYPRE", "MKLPDSS" ]:
            yield run, etype, solver
        if ESPRESOTest.fast: return

def run(etype, solver):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[7] = solver
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
