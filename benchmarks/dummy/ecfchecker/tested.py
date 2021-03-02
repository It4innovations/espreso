
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)

def teardown():
    os.remove(os.path.join(ESPRESOTest.path, "espreso.ecf.default"))

@istest
def printdefault():
    program = [ "{0}/build/ecfchecker".format(ESPRESOTest.root), "-d" ]
    output, error = ESPRESOTest.run_program(program)
    if error != "":
        ESPRESOTest.raise_error(error)
