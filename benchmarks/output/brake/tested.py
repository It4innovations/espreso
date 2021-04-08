
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "writers", "format", "writer" ]
    ESPRESOTest.store_results = True
    ESPRESOTest.external = True

def teardown():
    ESPRESOTest.clean("store")
    ESPRESOTest.clean("results")

@istest
def by():
    for processes in range(1, 32):
        for writers in [ 1, 12 ]:
            if writers <= processes:
                for writer in [ "MPI", "MPI_COLLECTIVE" ]:
                    for format in [ ("ENSIGHT", "store.case"), ("VTK_LEGACY", "store.*.vtk"), ("XDMF", "store.xmf") ]:
                        yield run, processes, writers, writer, format

def run(processes, writers, writer, format):
    ESPRESOTest.processes = processes
    ESPRESOTest.ecf = "store.ecf"
    ESPRESOTest.args[0] = writers
    ESPRESOTest.args[1] = writer
    ESPRESOTest.args[2] = format[0]
    ESPRESOTest.run()

    ESPRESOTest.processes = 16
    ESPRESOTest.ecf = "load.ecf"
    ESPRESOTest.args[0] = format[0]
    ESPRESOTest.args[1] = "store/last/{0}".format(format[1])
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
