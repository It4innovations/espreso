
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 2, 2, 2, 2, 20, 20, "COMPLEX" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "SQUARE4", "SQUARE8", "TRIANGLE3", "TRIANGLE6" ]:
        for system in [ "COMPLEX", "REAL" ]:
            yield run, etype, system

def run(etype, system):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[7] = system
    ESPRESOTest.run()

    if ESPRESOTest.create:
        ESPRESOTest.create_emr(".".join([etype, "emr"]))
    else:
        ESPRESOTest.compare_emr(".".join([etype, "emr"]))
