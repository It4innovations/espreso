
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.processes = 2

def teardown():
    ESPRESOTest.clean()

@istest
def defaults():
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
