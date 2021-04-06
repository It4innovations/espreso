
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "file", "decomposer", "readers", "loader" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for processes in range(2, 16):
        for decomposer in [ "PTSCOTCH", "HILBERT_CURVE" ]:
            for readers in [ 6, 15 ]:
                for loader in [ "MPI", "MPI_COLLECTIVE", "POSIX" ]:
                    yield run, processes, decomposer, readers, loader

def run(processes, decomposer, readers, loader):
    ESPRESOTest.processes = processes
    ESPRESOTest.args[0] = os.path.join(ESPRESOTest.path, "square.*.vtk")
    ESPRESOTest.args[1] = decomposer
    ESPRESOTest.args[2] = readers
    ESPRESOTest.args[3] = loader
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
