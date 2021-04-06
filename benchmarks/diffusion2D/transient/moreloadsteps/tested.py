
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ 2, 2, 3, 2, 10, 15, "solver" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for solver in [ "FETI", "HYPRE", "MKLPDSS" ]:
        yield run, solver

def run(solver):
    ESPRESOTest.args[6] = solver
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
    ESPRESOTest.report("espreso.time.xml")