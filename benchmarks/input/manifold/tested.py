
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
    path = os.path.join("/", "data", "espreso", "mesiotest", "manifold")
    ansys = os.path.join(path, "manifold.dat")
    ensightBinary = os.path.join(path, "manifold.case")
    ensightAscii = os.path.join(path, "manifoldAscii.case")
    vtk = os.path.join(path, "manifold.*.vtk")
    xdmf = os.path.join(path, "manifold.xmf")

    os.path.exists(ansys)
    os.path.exists(ensightBinary)
    os.path.exists(ensightAscii)
    os.path.exists(vtk)
    os.path.exists(xdmf)
    for p in range(8, 24, 3):
        for format, file in [ ("ANSYS_CDB", ansys), ("ENSIGHT", ensightBinary), ("ENSIGHT", ensightAscii), ("VTK_LEGACY", vtk), ("XDMF", xdmf) ]:
            for readers in [ 8, 19 ]:
                if readers < p:
                    yield run, file, format, p, readers, "MPI_COLLECTIVE"

def run(file, format, p, readers, loader):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = format
    ESPRESOTest.args[1] = file
    ESPRESOTest.args[2] = loader
    ESPRESOTest.args[3] = readers
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
