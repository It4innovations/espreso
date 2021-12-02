
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 2, 2, 1, 2, 2, 2, 4, 4, 8, "COMPLEX" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA8", "HEXA20", "TETRA4", "TETRA10", "PRISMA6", "PRISMA15", "PYRAMID5", "PYRAMID13" ]:
        for system in [ "COMPLEX", "REAL" ]:
            yield run, etype, system

def run(etype, system):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[10] = system
    ESPRESOTest.run()

    if ESPRESOTest.create:
        ESPRESOTest.create_emr(".".join([etype, "emr"]))
    else:
        ESPRESOTest.compare_emr(".".join([etype, "emr"]))
