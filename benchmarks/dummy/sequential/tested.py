
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "etype", 1, 1, 2, 2, 5, 5, "method" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for etype in [ "SQUARE4", "SQUARE8", "TRIANGLE3", "TRIANGLE6" ]:
        for method in [ "TOTAL_FETI", "HYBRID_FETI" ]:
            yield run, etype, method
            return

def run(etype, method):
    ESPRESOTest.args[0] = etype
    ESPRESOTest.args[7] = method
    ESPRESOTest.run()
    ESPRESOTest.report("espreso.time.xml")
