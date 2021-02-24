
import os
from nose.tools import istest

from estest import ESPRESOTest

def setup():
    ESPRESOTest.path = os.path.dirname(__file__)
    ESPRESOTest.args = [ "DIRICHLET", "ANALYTIC", "regularization_type", "conjugate_projector" ]

def teardown():
    ESPRESOTest.clean()

@istest
def by():
#     for regularization in [ "ANALYTIC", "ALGEBRAIC" ]:
    for rtype in [ "FIX_POINTS", "EIGEN_VECTORS" ]: # "WAVE_DIRECTIONS"
        for projection in [ "NONE", "CONJ_K" ]:
            yield run, rtype, projection

def run(rtype, projection):
#     ESPRESOTest.args[0] = preconditioner
#     ESPRESOTest.args[1] = regularization
    ESPRESOTest.args[2] = rtype
    ESPRESOTest.args[3] = projection
    ESPRESOTest.run()
    ESPRESOTest.compare_emr("espreso.emr")
