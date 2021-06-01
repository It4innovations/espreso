
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for processes in range(1, 16):
        yield run, processes

def run(processes):
    ESPRESOTest.processes = processes
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())

