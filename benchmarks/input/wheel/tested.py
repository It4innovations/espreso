
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
    path = os.path.join("/", "data", "espreso", "mesiotest", "wheel")
    ansys = os.path.join(path, "wheel.dat")
    ensightBinary = os.path.join(path, "wheel.case")
    ensightAscii = os.path.join(path, "wheelAscii.case")
    vtk = os.path.join(path, "wheel.*.vtk")
    xdmf = os.path.join(path, "wheel.xmf")

    os.path.exists(ansys)
    os.path.exists(ensightBinary)
    os.path.exists(ensightAscii)
    os.path.exists(vtk)
    os.path.exists(xdmf)
    for p in range(3, 24, 4):
        for format, file in [ ("ANSYS_CDB", ansys), ("ENSIGHT", ensightBinary), ("ENSIGHT", ensightAscii), ("VTK_LEGACY", vtk), ("XDMF", xdmf) ]:
            for readers in [ 3, 15 ]:
                if readers < p:
                    yield run, file, format, p, readers, "MPI"

def run(file, format, p, readers, loader):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = format
    ESPRESOTest.args[1] = file
    ESPRESOTest.args[2] = loader
    ESPRESOTest.args[3] = readers
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
