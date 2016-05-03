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

class Configuration:

    def __init__(self, name, parameters={}):
        self.name = name
        self.config = {
            "SUBDOMAINS": 8,
            "FIXPOINTS": 8,
            "CORNERS": 1,
            "VERTEX_CORNERS": 1,
            "EDGE_CORNERS": 1,

            "EPSILON": 1e-5,
            "ITERATIONS": 100,
            "FETI_METHOD": 0,
            "PRECONDITIONER": 1,
            "REGULARIZATION": 0,

            "REDUNDANT_LAGRANGE": 1,
            "USE_SCHUR_COMPLEMENT": 0,
            "KEEP_FACTORS": 1,
            "CGSOLVER": 0,
            "KSOLVER": 0,
            "KSOLVER_SP_iter_steps": 0,
            "KSOLVER_SP_iter_norm": 0,
            "F0SOLVER": 0,
            "N_MICS": 0,

            "FACE_CORNERS": 0,
            "AVERAGE_EDGES": 0,
            "AVERAGE_FACES": 0,

            "SAVE_MESH": 0,
            "SAVE_FIXPOINTS": 0,
            "SAVE_FACES": 0,
            "SAVE_EDGES": 0,
            "SAVE_CORNERS": 0,
            "SAVE_DIRICHLET": 0,
            "SAVE_AVERAGING": 0,
            "SAVE_RESULTS": 0,

            "VERBOSE_LEVEL": 0,
            "TESTING_LEVEL": 0,
            "MEASURE_LEVEL": 0,
            "PRINT_MATRICES": 0
        }

        for key in parameters:
            self.config[key] = parameters[key]

    def set(self, parameter, value):
        self.config[parameter] = value

    def save(self, path, suffix):
        file = open(os.path.join(ESPRESO_ROOT, EXAMPLES, path, "espreso.config.{0}.{1}".format(self.name, suffix)), 'w')
        for key, value in self.config.items():
            file.write("{0} = {1}\n".format(key, value))
        file.close()
        return "espreso.config.{0}.{1}".format(self.name, suffix)

class Espreso:

    def __init__(self, path, input, example, config):
        self.path = path
        self.input = input
        self.example = example
        self.config = config

    def run(self, processes, args, env=ENV):
        program = [ "mpirun", "-n", str(processes), os.path.join(ESPRESO_ROOT, "espreso")]
        program += [ "-p", self.example ]
        program += [ "-i", self.input ]
        program += [ "-c", self.config ]
        program += [ str(x) for x in args ]

        result = subprocess.Popen(program,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  cwd=os.path.join(EXAMPLES, self.path),
                                  env=dict(os.environ.items() + env.items()))

        return result.communicate()








