import os
from nose.tools import istest
from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__) # test will be executed from this directory
    # ESPRESOTest.args = [ "solver", "method" ] # test has 2 command line argument

def teardown():
    ESPRESOTest.clean()

@istest
def by():
    yield run

def run():
    ESPRESOTest.env["MKL_NUM_THREADS"]                   = "1"
    ESPRESOTest.env["OMP_NUM_THREADS"]                   = "1"
    ESPRESOTest.env["SOLVER_NUM_THREADS"]                = "1"
    ESPRESOTest.env["PAR_NUM_THREADS"]                   = "1"
    ESPRESOTest.env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("autoopt.emr")