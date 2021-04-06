
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "file", "decomposer", "readers" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for processes in range(1, 16):
        for decomposer in [ "PARMETIS", "PTSCOTCH" ]:
            for readers in [ 2, 3]:
                yield run, processes, decomposer, readers

def run(processes, decomposer, readers):
    ESPRESOTest.processes = processes
    ESPRESOTest.args[0] = os.path.join(ESPRESOTest.path, "cube.xmf")
    ESPRESOTest.args[1] = decomposer
    ESPRESOTest.args[2] = readers
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
