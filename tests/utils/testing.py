
import subprocess
import os
import sys
import re
import shutil

ESPRESO_UTILS = os.path.dirname(os.path.abspath(__file__))
ESPRESO_TESTS = os.path.dirname(ESPRESO_UTILS)
ESPRESO_ROOT = os.path.dirname(ESPRESO_TESTS)
EXAMPLES = os.path.join(ESPRESO_TESTS, "examples")

ENV = {
    "MKL_NUM_THREADS": "1",
    "CILK_NWORKERS": "1",

    "OMP_NUM_THREADS": "4",
    "SOLVER_NUM_THREADS": "4",
    "PAR_NUM_THREADS": "4",

    "PARDISOLICMESSAGE": "1",

    "OMP_PROC_BIND": "TRUE", # there is problem with linking CUDA
}


class TestCaseCreator:

    @staticmethod
    def create_test(testClass, function, name, *args, **kwargs):
        def test_method(self):
            function(self, *args, **kwargs)

        setattr(testClass, 'test_' + name, test_method)
        test_method.__name__ = 'test_' + name

    @staticmethod
    def iterate(function, *args):
        next = [ True for arg in args]
        iterators = [ Iterator(arg) for arg in args ]

        while reduce(lambda x, y: x or y, next):
            function(*[ it.get() for it in iterators ])
            for i in range(0, len(next)):
                next[i] = iterators[i].next()
                if next[i]:
                    break

    @staticmethod
    def select(*directories):
        argv = [ sys.argv[0] ]
        origin = []
        selection = []

        for directory in directories:
            origin.append([ directory ])
            selection.append([])
            for arg in sys.argv[1:]:
                if os.path.commonprefix([ directory, os.path.join(ESPRESO_ROOT, arg) ]) == directory:
                    selection[-1].append(os.path.join(ESPRESO_ROOT, arg))

        for arg in sys.argv[1:]:
            if arg.startswith("-") or arg.startswith("--"):
                argv.append(arg)

        if len(sys.argv) == len(argv):
            if len(origin) == 1:
                return origin[0]
            return tuple(origin)
        sys.argv = argv
        if len(selection) == 1:
            return selection[0]
        return tuple(selection)

    @staticmethod
    def gather(folder, ext, omit = "$^"):
        def skip(files):
            for file in files:
                if file.endswith(".skip"):
                    print "skip folder " + root
                    for line in open(os.path.join(root, file)):
                        print "--> " + line
                    return True
            return False

        # recurse tree a gather all examples
        omit = re.compile(omit)
        examples = []
        for root, subFolders, files in os.walk(folder):
            if "results" in root.split("/"):
                continue
            if not skip(files):
                for file in files:
                    if file.endswith(ext) and not omit.match(file):
                        examples.append(( root, file.replace(ext, '')))

        examples.sort()

        result = []
        for example in examples:
            name = os.path.relpath(example[0], folder).replace('/', '_')
            if example[1] != "espreso":
                name += "_" + example[1]
            result.append((name, example[0], example[1] + ext))
        return result

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

class EspresoError:

    program = []

    def __init__(self, error):
        raise Exception("{0}\nProgram: {1}".format(error, self.program))

class Espreso:

    def __init__(self, path=ESPRESO_ROOT, config={}):
        self.path = path
        self.mpirun = [ "mpirun" ]

        bashrc = os.path.join(os.path.expanduser("~"), ".bashrc")
        if os.path.isfile(bashrc):
            for line in open(bashrc, 'r'):
                if line.find("alias") != -1 and line.find("mpirun") != -1:
                    self.mpirun = map(lambda x: x.strip(" '\""), line.split("=")[1].split())

        if path != ESPRESO_ROOT:
            self.path = os.path.join(ESPRESO_TESTS, path)

        self.env = os.environ
        self.append_env("LD_LIBRARY_PATH", os.path.join(self.path, "libs/"))
        self.append_env("LD_LIBRARY_PATH", "/usr/local/cuda-7.5/lib64")
        for key, value in ENV.items():
            self.set_env(key, value)

        if path != ESPRESO_ROOT:
            self.clone()

    def append_env(self, key, value):
        self.env[key] = value + ":" + self.env[key]

    def set_env(self, key, value):
        self.env[key] = value


    def clone(self):
        self.remove()

        subprocess.call(["git", "clone", "git@code.it4i.cz:mec059/espreso.git", self.path, "-q"], cwd=ESPRESO_TESTS)
        subprocess.call(["mkdir", "libs"], cwd=self.path)

        origin = os.path.join(ESPRESO_ROOT, "libs/")
        target = os.path.join(self.path, "libs/")
        if os.path.isfile(origin + "libpardiso500-INTEL120-X86-64.so"):
            subprocess.call(["cp", origin + "libpardiso500-INTEL120-X86-64.so", target])
        if os.path.isfile(origin + "libifcore.a"):
            subprocess.call(["cp", origin + "libifcore.a", target])

    def remove(self):
        shutil.rmtree(self.path, ignore_errors=True)


    def waf(self, args=[], config={}):
        return self.run_program([os.path.join(self.path, "waf")], self.path, config, args)

    def install(self, config={}):
        def check(result, error, method):
            success = False
            for line in result.splitlines():
                if line.find("'" + method + "' finished successfully") != -1:
                    success = True

            if error != "":
                raise EspresoError(error)
            if success == False:
                raise EspresoError(result)

        result, error = self.waf(["configure"], config)
        check(result, error, "configure")
        result, error = self.waf(["install"])
        check(result, error, "install")

    def get_processes(self, ecf):
        def get_path(section):
            in_section = False
            for line in open(ecf, 'r'):
                if section in line.upper():
                    in_section = True
                if in_section and "PATH" in line.upper():
                    return line.split()[1].strip(";")

        def get_value(value):
            if value.startswith("[") and value.endswith("]"):
                return variables[value]
            return value

        variables = {}
        variable = False
        default_arg = False
        for line in open(ecf, 'r'):
            if "INPUT" in line.upper():
                input = line.split()[1].strip(";")

            if variable and "}" in line:
                variable = False
            if variable:
                variables["[" + line.split()[0].strip() + "]"] = line.split()[1].strip()
            if "VARIABLES" in line.upper():
                variable = True

            if default_arg and "}" in line:
                default_arg = False
            if default_arg and len(line.strip()):
                variables["[ARG" + line.split()[0].strip() + "]"] = line.split()[1].strip(";")
            if "DEFAULT_ARGS" in line.upper():
                default_arg = True

        if input.upper() == "GENERATOR":
            in_generator = 0
            in_tower = 0
            tower = False
            sphere = False
            multiplier = 1
            proc_index = 0
            parrents = 0
            for line in open(ecf, 'r'):
                if in_generator and "{" in line:
                    in_generator += 1
                if in_tower and "{" in line:
                    in_tower += 1
                if "GENERATOR" in line.upper():
                    in_generator = 1
                if in_generator and "GRIDS" in line.upper():
                    in_tower = 1
                if in_generator and "}" in line:
                    in_generator -= 1
                if in_tower and "}" in line:
                    in_tower -= 1
                if in_generator and "GRID" in line.upper():
                    procs = [ (1, 1, 1) ]
                if in_generator and "GRID_TOWER" in line.upper():
                    procs = [ ]
                    tower = True
                if in_generator and "SPHERE" in line.upper():
                    procs = [ (1, 1, 1) ]
                    sphere = True
                    multiplier = 6

                if in_generator and in_tower and len(line.strip()) and line.split()[0].isdigit():
                    proc_index = int(line.split()[0])
                    if proc_index >= len(procs):
                        if proc_index != len(procs):
                            EspresoError("Not supported format of *ecf test file. Default GRIDS settings.")
                        procs.append((1, 1, 1))
                if in_generator and "CLUSTERS_X" in line.upper():
                    procs[proc_index] = (int(get_value(line.split()[1].strip(";"))), procs[proc_index][1], procs[proc_index][2])
                if in_generator and "CLUSTERS_Y" in line.upper():
                    procs[proc_index] = (procs[proc_index][0], int(get_value(line.split()[1].strip(";"))), procs[proc_index][2])
                if in_generator and "CLUSTERS_Z" in line.upper():
                    procs[proc_index] = (procs[proc_index][0], procs[proc_index][1], int(get_value(line.split()[1].strip(";"))))
                if sphere and "CLUSTERS" in line.upper():
                    procs[proc_index] = (int(get_value(line.split()[1].strip(";"))), procs[proc_index][1], procs[proc_index][2])
                if sphere and "LAYERS" in line.upper():
                    procs[proc_index] = (procs[proc_index][0], int(get_value(line.split()[1].strip(";"))), procs[proc_index][2])

            return multiplier * sum(i[0] * i[1] * i[2] for i in procs)

        if input.upper() == "OPENFOAM":
            return len(glob.glob(get_path("OPENFOAM") + "/processor*"))

        if input.upper() == "ANSYS":
            return 1

        if input.upper() == "ESDATA":
            return len(filter(os.path.isdir, glob.glob(get_path("ESDATA") + "/*")))

    def compare_monitors(self, emr1, emr2):
        def compare(v1, v2):
            if v1 != v2 and abs(float(v1) - float(v2)) > 1e-4 and abs((float(v1) - float(v2)) / float(v1)) > 1e-3:
                raise EspresoError(
                    "various monitored results" +
                    "[" + table1[0][column] + "::" + table1[1][column] + "::" + table1[2][column] + "]" +
                    "[step " + table1[row][0] + "][substep " + table1[row][1] + "]" +
                    " -> " + v1 + " != " + v2
                )

        def create_table(emr):
            return [ [ value.strip() for value in line.split(";") ] for line in open(emr, "r").readlines() if len(line.strip()) ]

        table1 = create_table(emr1)
        table2 = create_table(emr2)

        if len(table1[0]) != len(table2[1]):
            raise EspresoError("various monitored properties")

        for row, (row1, row2) in enumerate(zip(table1, table2)):
            for column, (value1, value2) in enumerate(zip(row1, row2)):
                compare(value1, value2)

    def run_program(self, program, cwd="", config={}, args=[]):
        config["ENV::REMOVE_OLD_RESULTS"] = "1"
        program += [ str(x) for x in args ]
        for key, value in config.items():
            program += [ "--{0}={1}".format(key, value) ]

        if cwd:
            cwd = os.path.join(ESPRESO_ROOT, cwd)
        EspresoError.program = program
        result = subprocess.Popen(program,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  cwd=cwd or self.path,
                                  env=self.env)

        return result.communicate()


    def run(self, processes, *args, **kwargs):
        program = self.mpirun + [ "-n", str(processes),"--map-by", "core", os.path.join(self.path, "espreso")]

        output, error = self.run_program(program, *args, **kwargs)
        if error != "":
            raise EspresoError(error)
        if output != "":
            raise EspresoError(output)

    def valgrind(self, processes, *args, **kwargs):
        program = self.mpirun + [ "-n", str(processes), "valgrind", "-q", "--leak-check=full", "--suppressions={0}/espreso.supp".format(self.path), os.path.join(self.path, "espreso")]

        output, error = self.run_program(program, *args, **kwargs)
        if error != "":
            skip = False
            warningless = ""
            for line in error.split("\n"):
                tokens = line.split(" ")
                if len(tokens) == 1:
                    continue
                if tokens[1] == "Warning:":
                    skip = True
                    continue
                elif tokens[1] == "" and skip:
                    continue
                else:
                    skip = False
                warningless += line + "\n"
            if warningless:
                raise EspresoError("\n" + warningless)
        if output != "":
            raise EspresoError(output)

    def decompose(self, *args, **kwargs):
        program = [ os.path.join(self.path, "decomposer") ]

        output, error = self.run_program(program, *args, **kwargs)
        if error != "":
            raise EspresoError(error)
        return output

    def output(self, processes, *args, **kwargs):
        program = self.mpirun + [ "-n", str(processes), os.path.join(self.path, "espreso")]

        output, error = self.run_program(program, *args, **kwargs)
        if error != "":
            raise EspresoError(error)

        return output

    def fail(self, processes, *args, **kwargs):
        program = self.mpirun + [ "-n", str(processes), os.path.join(self.path, "espreso")]

        output, error = self.run_program(program, *args, **kwargs)
        if error == "":
            raise EspresoError("Expected fail, but run was correct.")






