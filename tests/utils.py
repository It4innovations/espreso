import subprocess
import os

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ESPRESO_ROOT = os.path.dirname(ESPRESO_TESTS)
EXAMPLES = os.path.join(ESPRESO_TESTS, "examples")

ENV = {
    "CILK_NWORKERS": "2",
    "SOLVER_NUM_THREADS": "2",
    "MKL_NUM_THREADS": "2",
    "PAR_NUM_THREADS": "2",
    "OMP_NUM_THREADS": "2"
}

if not os.path.isdir(EXAMPLES):
    os.mkdir(EXAMPLES)

class Espreso:

    def __init__(self, example):
        self.example = EXAMPLES + "/" + example

    def run(self, processes, args, env=ENV):
        program = [ "mpirun", "-n", str(processes), os.path.join(ESPRESO_ROOT, "espreso")]
        program += [ str(x) for x in args]

        result = subprocess.Popen(program,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  cwd=self.example,
                                  env=dict(os.environ.items() + env.items()))

        return result.communicate()








