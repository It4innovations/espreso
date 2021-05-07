
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 1, 1, 1, 3, 2, 4, 3, 5, 2 ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA8", "HEXA20", "TETRA4", "TETRA10", "PRISMA6", "PRISMA15", "PYRAMID5", "PYRAMID13" ]:
        yield run, etype
        if ESPRESOTest.fast: return

def run(etype):
    ESPRESOTest.processes = 1
    ESPRESOTest.args[0] = etype
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))