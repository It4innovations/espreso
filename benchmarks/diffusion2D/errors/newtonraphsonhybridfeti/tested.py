
import os
from nose.tools import istest
from unittest.case import skip

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.getcwd()
    ESPRESOTest.args = []

def teardown():
    ESPRESOTest.clean()

@istest
@skip("HYBRID FETI makes some illness")
def defaults():
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
    ESPRESOTest.report("espreso.time.xml")
