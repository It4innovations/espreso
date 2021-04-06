
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 1, 1, 2, 1, 1, 2, 2, 2, 10, "solver" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    yield run, "HEXA8", "MKLPDSS"
    yield run, "HEXA20", "MKLPDSS"
    yield run, "HEXA8", "FETI"

def run(etype, solver):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[10] = solver
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
