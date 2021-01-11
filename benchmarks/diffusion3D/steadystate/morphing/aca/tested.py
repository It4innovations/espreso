
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "file" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for p in range(1, 32):
        yield run, p
        if ESPRESOTest.fast: return

def run(p):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = os.path.join(ESPRESOTest.path, "espreso.case")
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")