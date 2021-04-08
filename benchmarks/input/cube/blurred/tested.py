
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
    for processes in range(1, 16):
        yield run, processes

def run(processes):
    ESPRESOTest.processes = processes
    ESPRESOTest.args[0] = os.path.join(ESPRESOTest.path, "blurred.case")
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
