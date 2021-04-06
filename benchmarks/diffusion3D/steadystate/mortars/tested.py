
import os
from nose.tools import istest

from estest import ESPRESOTest
from unittest.case import skip

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.processes = 2
    ESPRESOTest.args = []

def teardown():
    ESPRESOTest.clean()

@istest
@skip("to be improved soon")
def defaults():
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
    ESPRESOTest.report("espreso.time.xml")
