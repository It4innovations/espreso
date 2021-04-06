
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 1, 1, 4, 2, 2, 2, 2, 2, 5 ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "HEXA8", "HEXA20" ]:
        yield run, etype
        if ESPRESOTest.fast: return

def run(etype):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.run()
    ESPRESOTest.compare_emr(".".join([etype, "emr"]))
    ESPRESOTest.report("espreso.time.xml")
