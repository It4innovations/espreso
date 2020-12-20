
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "loader" ]
    ESPRESOTest.external = True

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for loader in [ "POSIX", "MPI", "MPI_COLLECTIVE" ]:
        yield run, 1, loader
        yield run, 2, loader

def run(p, loader):
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = loader
    ESPRESOTest.compare_mesh("espreso.log", ESPRESOTest.run())
