
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.env["OMP_NUM_THREADS"] = "4"
    ESPRESOTest.env["SOLVER_NUM_THREADS"] = "4"
    ESPRESOTest.env["PAR_NUM_THREADS"] = "4"
    ESPRESOTest.args = [ 1, 20, "TRUE", 1.5, 0]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    for p in [1, 2, 4]:
        for sc_and_spds in ["TRUE" , "FALSE"]:
            for fits_on_gpu in ["ALL", "PARTIAL", "NONE"]:
		yield run, p, sc_and_spds, fits_on_gpu

def run(p, sc_and_spds, fits_on_gpu):
    
    ESPRESOTest.processes = p
    ESPRESOTest.args[0] = p 
    ESPRESOTest.args[1] = 20/p 
    ESPRESOTest.args[2] = sc_and_spds
    if fits_on_gpu == "PARTIAL":
	if p == 1:
            ESPRESOTest.args[3] = 8
	if p == 2:
            ESPRESOTest.args[3] = 25
        if p == 4:
            ESPRESOTest.args[3] = 256
    	ESPRESOTest.args[4] = 256
    if fits_on_gpu == "NONE":
    	ESPRESOTest.args[3] = 2000 
    	ESPRESOTest.args[4] = 256
    ESPRESOTest.run()	
    ESPRESOTest.compare_emr("espreso.emr")
