
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.store_results = True

def teardown():
    ESPRESOTest.clean("store")
    ESPRESOTest.clean("results")

@istest
def by():
    for p in range(1, 32):
        yield run, p

def run(p):
    ESPRESOTest.processes = 1
    ESPRESOTest.ecf = "store.ecf"
    ESPRESOTest.run()

    ESPRESOTest.processes = p
    ESPRESOTest.ecf = "load.ecf"
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
