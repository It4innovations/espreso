
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "format", "file", "loader", "readers" ]
    ESPRESOTest.external = True

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    path = os.path.join("/", "data", "espreso", "mesiotest", "cooler", "linear")
    ansys = os.path.join(path, "cooler.dat")
    ensightBinary = os.path.join(path, "cooler.case")
    ensightAscii = os.path.join(path, "coolerAscii.case")
    vtk = os.path.join(path, "cooler.*.vtk")
    xdmf = os.path.join(path, "cooler.xmf")

    os.path.exists(ansys)
    os.path.exists(ensightBinary)
    os.path.exists(ensightAscii)
    os.path.exists(vtk)
    os.path.exists(xdmf)
    for p in range(1, 32, 4):
        for format, file in [ ("ANSYS_CDB", ansys), ("ENSIGHT", ensightBinary), ("ENSIGHT", ensightAscii), ("VTK_LEGACY", vtk), ("XDMF", xdmf) ]:
            for readers in [ 2, 7, 12, 23 ]:
                if readers < p:
                    yield run, file, format, p, readers, [ "POSIX", "MPI", "MPI_COLLECTIVE" ][int(((p + 1) / 3) % 3)]

def run(file, format, p, readers, loader):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = format
    ESPRESOTest.args[1] = file
    ESPRESOTest.args[2] = loader
    ESPRESOTest.args[3] = readers
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
