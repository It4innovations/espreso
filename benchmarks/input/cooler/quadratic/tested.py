
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "format", "file", "loader", "readers" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    path = os.path.join("/", "data", "espreso", "mesiotest", "cooler", "quadratic")
    ansys = os.path.join(path, "cooler.dat")
    ensight = os.path.join(path, "cooler.case")
    vtk = os.path.join(path, "cooler.*.vtk")
    xdmf = os.path.join(path, "cooler.xmf")

    os.path.exists(ansys)
    os.path.exists(ensight)
    os.path.exists(vtk)
    os.path.exists(xdmf)
    for p in range(4, 32, 6):
        for format, file in [ ("ANSYS_CDB", ansys), ("ENSIGHT", ensight), ("VTK_LEGACY", vtk), ("XDMF", xdmf) ]:
                for readers in [ 4, 13, 29 ]:
                    if readers < p:
                        yield run, file, format, p, readers, "MPI_COLLECTIVE"
                    if readers == p:
                        for loader in [ "MPI_COLLECTIVE", "MPI", "POSIX" ]:
                            yield run, file, format, p, readers, loader

def run(file, format, p, readers, loader):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = format
    ESPRESOTest.args[1] = file
    ESPRESOTest.args[2] = loader
    ESPRESOTest.args[3] = readers
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
