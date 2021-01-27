
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 1, 1, 1, 3, 2, 4, 16, 16, 32 ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA8" ]:
        for p in range(1, 32):
            yield run, p, etype
            if ESPRESOTest.fast: return

def run(p, etype):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[3] = p
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))