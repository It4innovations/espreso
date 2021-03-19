
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ 1, 8 ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for p in [1, 2, 4, 8]:
	for t in ["1", "2", "4"]:
		yield run, p, t

def run(p, t):
    ESPRESOTest.env["OMP_NUM_THREADS"] = t
    ESPRESOTest.env["SOLVER_NUM_THREADS"] = t
    ESPRESOTest.env["PAR_NUM_THREADS"] = t
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = p 
    ESPRESOTest.args[1] = 8/p 
    ESPRESOTest.run()	
    ESPRESOTest.compare_emr("espreso.emr")
