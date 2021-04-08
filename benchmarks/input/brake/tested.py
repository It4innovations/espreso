
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
    path = os.path.join("/", "data", "espreso", "mesiotest", "brake")
    ansys = os.path.join(path, "brake.dat")
    ensight = os.path.join(path, "brake.case")
    vtk = os.path.join(path, "brake.*.vtk")
    xdmf = os.path.join(path, "brake.xmf")

    os.path.exists(ansys)
    os.path.exists(ensight)
    os.path.exists(vtk)
    os.path.exists(xdmf)
    for p in range(1, 32, 3):
        for format, file in [ ("ANSYS_CDB", ansys), ("ENSIGHT", ensight), ("VTK_LEGACY", vtk), ("XDMF", xdmf) ]:
            for readers in [ 1, 7, 19, 29 ]:
                if readers < p:
                    yield run, file, format, p, readers, [ "POSIX", "MPI", "MPI_COLLECTIVE" ][int((p / 3) % 3)]

def run(file, format, p, readers, loader):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = format
    ESPRESOTest.args[1] = file
    ESPRESOTest.args[2] = loader
    ESPRESOTest.args[3] = readers
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
