
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 2, 2, 1, 2, 2, 2, 4, 4, 8 ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA8", "HEXA20", "TETRA4", "TETRA10", "PRISMA6", "PRISMA15", "PYRAMID5", "PYRAMID13" ]:
        yield run, etype

def run(etype):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("{}.emr".format(etype))
