
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
    path = os.path.join("/", "data", "espreso", "mesiotest", "drillbit")
    ansys = os.path.join(path, "drillbit.dat")
    ensight = os.path.join(path, "drillbit.case")
    vtk = os.path.join(path, "drillbit.*.vtk")
    xdmf = os.path.join(path, "drillbit.xmf")

    os.path.exists(ansys)
    os.path.exists(ensight)
    os.path.exists(vtk)
    os.path.exists(xdmf)
    for p in range(2, 16, 2):
        for format, file in [ ("ANSYS_CDB", ansys), ("ENSIGHT", ensight), ("VTK_LEGACY", vtk), ("XDMF", xdmf) ]:
            for readers in [ 1, 3, 13 ]:
                if readers < p:
                    yield run, file, format, p, readers, "POSIX"

def run(file, format, p, readers, loader):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = format
    ESPRESOTest.args[1] = file
    ESPRESOTest.args[2] = loader
    ESPRESOTest.args[3] = readers
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
