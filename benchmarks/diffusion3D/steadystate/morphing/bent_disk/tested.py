
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 1, 1, 1, 3, 4, 2, 16, 16, 16 ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA8" ]:
        for p in [[1, 1, 1], [2, 1, 1], [2, 2, 1], [2, 2, 2], [4, 2, 2], [4, 4, 2]]:
            dim = 16;
            yield run, p[0], p[1], p[2], dim, etype
            if ESPRESOTest.fast: return

def run(p1, p2, p3, dim, etype):
    ESPRESOTest.processes = p1*p2*p3
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[1] = p1
    ESPRESOTest.args[2] = p2
    ESPRESOTest.args[3] = p3
    ESPRESOTest.args[7] = dim/p1
    ESPRESOTest.args[8] = dim/p2
    ESPRESOTest.args[9] = dim/p3
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))