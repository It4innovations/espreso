
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "writer" ]
    ESPRESOTest.store_results = True
    ESPRESOTest.external = True

def teardown():
    ESPRESOTest.clean("store")
    ESPRESOTest.clean("results")

@istest
def by():
    for processes in [1]:
        for writer in [ "MPI", "MPI_COLLECTIVE" ]:
                yield run, processes, writer

def run(processes, writer):
    ESPRESOTest.processes = processes
    ESPRESOTest.ecf = "store.ecf"
    ESPRESOTest.args[0] = writer
    ESPRESOTest.run()

    ESPRESOTest.processes = 16
    ESPRESOTest.ecf = "load.ecf"
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
