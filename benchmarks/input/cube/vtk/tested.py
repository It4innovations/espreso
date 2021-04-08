
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "file", "decomposer", "loader" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for processes in range(1, 16):
        for decomposer in [ "PARMETIS" ]:
            for loader in [ "MPI", "MPI_COLLECTIVE", "POSIX" ]:
                yield run, processes, decomposer, loader

def run(processes, decomposer, loader):
    ESPRESOTest.processes = processes
    ESPRESOTest.args[0] = os.path.join(ESPRESOTest.path, "cube.*.vtk")
    ESPRESOTest.args[1] = decomposer
    ESPRESOTest.args[3] = loader
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
