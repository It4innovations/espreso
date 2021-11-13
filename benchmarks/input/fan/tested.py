
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
    path = os.path.join("/", "data", "espreso", "mesiotest", "fan")
    ansys = os.path.join(path, "fan.dat")
    ensightAscii = os.path.join(path, "fanAscii.case")
    ensightBinary = os.path.join(path, "fan.case")
    vtk = os.path.join(path, "fan.*.vtk")
    xdmf = os.path.join(path, "fan.xmf")

    os.path.exists(ansys)
    os.path.exists(ensightAscii)
    os.path.exists(ensightBinary)
    os.path.exists(vtk)
    os.path.exists(xdmf)
    for p in range(1, 32):
        for format, file in [ ("ANSYS_CDB", ansys), ("ENSIGHT", ensightBinary), ("ENSIGHT", ensightAscii), ("VTK_LEGACY", vtk), ("XDMF", xdmf) ]:
            for readers in [ 2, 19 ]:
                if readers < p:
                    yield run, file, format, p, readers, "MPI"

def run(file, format, p, readers, loader):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = format
    ESPRESOTest.args[1] = file
    ESPRESOTest.args[2] = loader
    ESPRESOTest.args[3] = readers
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
