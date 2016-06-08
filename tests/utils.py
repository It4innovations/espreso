import subprocess
import os

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ESPRESO_ROOT = os.path.dirname(ESPRESO_TESTS)
EXAMPLES = os.path.join(ESPRESO_TESTS, "examples")

ENV = {
    "MKL_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",

    "CILK_NWORKERS": "2",
    "SOLVER_NUM_THREADS": "2",
    "PAR_NUM_THREADS": "2"
}


class TestCaseCreator:

    @staticmethod
    def create_test(testClass, function, name, *args, **kwargs):
        def test_method(self):
            function(self, *args, **kwargs)

        setattr(testClass, 'test_' + name, test_method)
        test_method.__name__ = 'test_' + name

    @staticmethod
    def iterate(function, *args, **kwargs):
        next = [ True for arg in args]

        while reduce(lambda x, y: x or y, next):
            function(*args, **kwargs)
            for i in range(0, len(next)):
                next[i] = args[i].next()
                if next[i]:
                    break


class Iterator:

    def __init__(self, items):
        self.items = items
        self.keys = items.keys()
        self.pointers = [ 0 for i in items ]

    def next(self):
        self.pointers[-1] += 1
        for i in xrange(1, len(self.pointers)):
            if self.pointers[-i] == len(self.items[self.keys[-i]]):
                self.pointers[-i] = 0
                self.pointers[-i - 1] += 1

        if self.pointers[0] == len(self.items[self.keys[0]]):
            self.pointers[0] = 0
            return False

        return True

    def __getitem__(self, i):
        return self.get()[i]

    def values(self):
        return self.get().values()

    def get(self):
        result = {}
        for i in xrange(0, len(self.pointers)):
            result[self.keys[i]] = self.items[self.keys[i]][self.pointers[i]]
        return result

class Espreso:

    def __init__(self, path):
        self.path = path
        self.root = ESPRESO_ROOT

    def run_program(self, program, args, config, env=ENV):
        program += [ str(x) for x in args ]
        for key, value in config.items():
            program += [ "--{0}={1}".format(key, value) ]

        env=dict(os.environ.items() + env.items())
        env["LD_LIBRARY_PATH"] = os.path.join(self.root, "libs/") + ":" + env["LD_LIBRARY_PATH"]
        result = subprocess.Popen(program,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  cwd=os.path.join(EXAMPLES, self.path), env=env)

        return result.communicate()

    def configure(self, config):
        self.root = os.path.join(EXAMPLES, self.path)
        program = [ os.path.join(self.root, "waf"), "configure"]
        return self.run_program(program, [], config)

    def build(self):
        program = [ os.path.join(self.root, "waf"), "install"]
        return self.run_program(program, [], {})

    def run(self, processes, example, input, args, config, env=ENV):
        program = [ "mpirun", "-n", str(processes), os.path.join(self.root, "espreso")]
        program += [ "-p", example ]
        program += [ "-i", input ]
        return self.run_program(program, args, config, env)








