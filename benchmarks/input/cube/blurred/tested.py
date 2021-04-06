
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "file", "readers" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for processes in range(1, 16):
        for readers in [ 1, 3, 15]:
            if readers <= processes:
                yield run, processes, readers

def run(processes, readers):
    ESPRESOTest.processes = processes
    ESPRESOTest.args[0] = os.path.join(ESPRESOTest.path, "blurred.case")
    ESPRESOTest.args[1] = readers
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
