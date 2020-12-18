
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for p in range(1, 4):
        yield run, p

def run(p):
    ESPRESOTest.processes = p
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
