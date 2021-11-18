
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "format", "file" ]
    ESPRESOTest.external = True

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    # VTK is not tested since rouding errors in ASCII format
    path = os.path.join("/", "data", "espreso", "mesiotest", "fan_part")
    ansys = os.path.join(path, "fan_part.dat")
    ensight = os.path.join(path, "fan_part.case")
    xdmf = os.path.join(path, "fan_part.xmf")

    os.path.exists(ansys)
    os.path.exists(ensightBinary)
    os.path.exists(ensightAscii)
    os.path.exists(xdmf)
    for p in range(8, 32):
        for format, file in [ ("ANSYS_CDB", ansys), ("ENSIGHT", ensight), ("XDMF", xdmf) ]:
            yield run, file, format, p

def run(file, format, p):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = format
    ESPRESOTest.args[1] = file
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
