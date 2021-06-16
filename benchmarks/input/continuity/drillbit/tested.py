
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.external = True

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for p in range(10, 33):
        yield run, p

def run(p):
    ESPRESOTest.processes = p
    ESPRESOTest.check_continuity(ESPRESOTest.run())
