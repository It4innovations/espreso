
import os
from nose.tools import istest

from estest import ESPRESOTest
from unittest.case import skip

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 1, 1, 2, 2, 2, 4, 4, 4, ]

def teardown():
    ESPRESOTest.clean()

@istest
@skip("Invalid sphere generator")
def by():
    for etype in [ "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ]:
        yield run, etype

def run(etype):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
