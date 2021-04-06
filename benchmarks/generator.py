
import os
import sys
import shutil
import time
import math
from compiler.syntax import check

PATH = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(PATH)
RUNDIR = os.getcwd()

class ECFParser:

    def tokenize(self, file):
        tokens = []
        comment = False
        for line in open(file, "r"):
            for token in line.replace(";", " ;").split():
                if not comment and (token == "#" or token == "//"):
                    break
                if token == "/*":
                    comment = True
                if not comment:
                    tokens.append(token)
                if token == "*/":
                    comment = False

        return tokens

    def get_object(self, tokens, name):
        object = []
        brackets = 0
        is_in = False
        for token in tokens:
            if token == name:
                is_in = True
            if is_in and token == "}":
                brackets -= 1
                if brackets == 0:
                    return object
            if brackets > 0:
                object.append(token)
            if is_in and token == "{":
                brackets += 1

    def get_ecf_tree(self, tokens):
        tree = {}
        stack = [tree]
        values = []
        for token in tokens:
            values.append(token)

            if values[-1] == "}":
                stack.pop()
                values = []
            if len(values) > 1 and values[-1] == "{":
                stack[-1][values[0].upper()] = {}
                stack.append(stack[-1][values[0].upper()])
                values = []
            if len(values) > 1 and values[-1] == ";":
                stack[-1][values[0].upper()] = " ".join(values[1:-1])
                values = []

        self.check_ecf(tree)
        return tree

    def print_ecf(self, ecf, indent):
        def spaces(indent):
            return "".join([" " for i in range(0, indent)])

        for parameter in ecf:
            if isinstance(ecf[parameter], dict):
                print spaces(indent) + parameter + " {"
                self.print_ecf(ecf[parameter], indent + 2)
                print spaces(indent) + "}"
            else:
                print spaces(indent) + parameter + " = " + ecf[parameter]

    def check_ecf(self, ecf):
        def check_parameter(object, name, error):
            if not name in object or not len(object[name]):
                print error
                exit()

        check_parameter(ecf, "OUTPUT", "Set PYTHON_TEST_GENERATOR::OUTPUT")
        check_parameter(ecf, "LEVELS", "Set at least one LEVEL")
        check_parameter(ecf, "ENV", "Set environment setting")
        check_parameter(ecf, "RUN", "Set run command")
        check_parameter(ecf, "EXE", "Set execute command")
        check_parameter(ecf, "LOG", "Set log directory")
        check_parameter(ecf, "GATHER_LEVEL", "Set the level for that the various execute files will be generated")

        if "TABLES" not in ecf:
            return

        levels = len(ecf["LEVELS"])
        for table in ecf["TABLES"]:
            check_parameter(ecf["TABLES"][table], "ROWS", "Set ROWS for table '{0}'".format(table))
            check_parameter(ecf["TABLES"][table], "COLUMNS", "Set COLUMNS for table '{0}'".format(table))
            check_parameter(ecf["TABLES"][table], "VALUES", "Set values for table '{0}'".format(table))
            for level in range(1, levels + 1):
                check_parameter(ecf["TABLES"][table], "L{0}".format(level), "Each table has to set values for all levels")

            for param in ecf["TABLES"][table]:
                if param not in ["ROWS", "COLUMNS", "VALUES"] + ["L{0}".format(i) for i in range(1, levels + 1)]:
                    print "Table {0} constains invalid parameter '{1}'".format(table, param)

class ESPRESOTestGenerator:

    def __init__(self):
        self.ROOT = ROOT
        self.TIME = time.strftime("%H-%M-%S")
        self.DATE = time.strftime("%Y-%m-%d")

    def create_dir(self, dir):
        if not os.path.isdir(dir):
            os.makedirs(dir)
        return dir

    def create_file(self, path, levels, extension):
        return open(os.path.join(path, ".".join([name for value, name in levels]) + extension), "w")

    def evaluate(self, expr, variables, path):
        mathfnc = dict((name, fnc) for name, fnc in math.__dict__.iteritems() if callable(fnc))
        try:
            result = eval(expr, variables, mathfnc)
            del variables["__builtins__"]
            return result
        except:
            print "Invalid expression 'PYTHON_TEST_GENERATOR::" + path + "'."
            exit()

    def generate_level(self, ecf, run, log, stats, levels=[]):
        run_files = []
        exe_files = []
        for value, name in self.evaluate(ecf["LEVELS"][str(len(levels) + 1)], {}, "PYTHON_TEST_GENERATOR::LEVELS::" + str(len(levels) + 1)):
            subpath = self.create_dir(os.path.join(run, name))
            levels.append((value, name))

            if str(len(levels) + 1) in ecf["LEVELS"]:
                run_file, executables = self.generate_level(ecf, subpath, log, stats, levels)
                run_files.append(run_file)
                exe_files += executables
            else:
                run_files.append(self.generate_example(ecf, subpath, log, levels))
                if len(levels) == int(ecf["GATHER_LEVEL"]):
                    exe_files.append((run_files[-1], [item for item in levels]))
            levels.pop()

        run_file = self.gather_run_files(ecf, run_files, run, levels)
        if len(levels) == int(ecf["GATHER_LEVEL"]):
            exe_files.append((run_file, [item for item in levels]))

        if len(levels) == 0:
            self.generate_execute_file(ecf, exe_files)
        return run_file, exe_files

    def generate_example(self, ecf, run, log, levels):
        variables = {
            "ROOT": self.ROOT,
            "TIME": self.TIME,
            "DATE": self.DATE
        }

        for idx, (value, name) in enumerate(levels):
            variables["L" + str(idx + 1)] = value

        if "VARIABLES" in ecf:
            for variable, expression in ecf["VARIABLES"].iteritems():
                variables[variable] = self.evaluate(expression, variables, "VARIABLES::" + variable)

        n = 1
        if "MEASURE_REPETITION" in ecf:
            n = int(ecf["MEASURE_REPETITION"])

        repetitions = []
        for i in range(0, n):
            example = self.create_file(run, levels, ".{0}.ecf".format(i))
            repetitions.append(example)
            example.write("# Automatically generated test\n\n")

            args = {}
            if "ARGS" in ecf:
                args = ecf["ARGS"]
            example.write("DEFAULT_ARGS { # automatically generated ARGUMENTS\n")
            for arg, value in sorted(args.iteritems()):
                value = self.evaluate(value, variables, "ARGS::" + arg)
                example.write("  {0} {1};\n".format(arg, value))
            example.write("} # end of generated ARGUMENTS\n\n")

            for line in open(ecf["origin"], "r"):
                example.write(line)

        run_file = self.create_file(run, levels, ".run.sh")
        os.chmod(run_file.name, 0o777)
        run_file.write("#!/bin/bash\n\n")
        run_file.write("{0}\n".format(self.evaluate(ecf["ENV"], variables, "ENV")))
        run_file.write("cd {0}\n\n".format(run))

        if "WARMUP" in ecf:
            run_file.write("# WARM UP\n")
            run_file.write("{0} -c {1}\n\n".format(self.evaluate(ecf["WARMUP"], variables, "WARMUP"), example.name))

        for example in repetitions:
            run_file.write("{0} -c {1}\n".format(self.evaluate(ecf["RUN"], variables, "RUN"), example.name))
            run_file.write("cp {0}/*.log {1}\n\n".format(self.evaluate(ecf["LOG"], variables, "LOG"), log))

        if "POST" in ecf:
            run_file.write("# POST\n")
            run_file.write("{0}\n\n".format(self.evaluate(ecf["POST"], variables, "POST"), example.name))

        return run_file.name

    def gather_run_files(self, ecf, run_files, path, levels):
        if len(levels):
            run_file = self.create_file(path, levels, ".run.sh")
        else:
            run_file = self.create_file(path, levels, "{0}.run.sh".format(ecf["OUTPUT"]))
        os.chmod(run_file.name, 0o777)
        for file in run_files:
            run_file.write("{0}\n".format(file))

        return run_file.name

    def generate_execute_file(self, ecf, run_files):
        def check_variables(expr, denied):
            passed = True
            for variable in [item for item in denied]:
                if variable in expr:
                    denied.append(variable)
                    passed = False
            return passed

        variables = {
            "ROOT": self.ROOT,
            "TIME": self.TIME,
            "DATE": self.DATE
        }

        denied = []
        for idx, value in enumerate(ecf["LEVELS"].iteritems()):
            if idx >= int(ecf["GATHER_LEVEL"]):
                denied.append("L" + str(idx + 1))

        execute_file = self.create_file(os.path.join(ROOT, "generatedtests", ecf["OUTPUT"]), [], "execute.sh")
        os.chmod(execute_file.name, 0o777)
        for file, levels in run_files:
            for idx, (value, name) in enumerate(levels):
                variables["L" + str(idx + 1)] = value
            if "VARIABLES" in ecf:
                for variable, expression in ecf["VARIABLES"].iteritems():
                    if check_variables(expression, denied):
                        variables[variable] = self.evaluate(expression, variables, "VARIABLES::" + variable)

            if check_variables(ecf["EXE"], denied):
                execute_file.write("{0} {1}\n".format(self.evaluate(ecf["EXE"], variables, "EXE"), file))
            else:
                print "Invalid GATHER_LEVEL: EXE command references variable that vary in sublevels."
                exit()

    def generate_evaluate_file(self, ecf, log, stats):
        if "TABLES" not in ecf:
            return

        variables = {
            "ROOT": self.ROOT,
            "TIME": self.TIME,
            "DATE": self.DATE
        }

        levels = len(ecf["LEVELS"])
        for level in range(1, levels + 1):
            values = []
            for (value, name) in self.evaluate(ecf["LEVELS"][str(level)], {}, ""): # error is irelevant, was evaluated before
                values.append("'{0}'".format(name))
            variables["L{0}".format(level)] = "[{0}]".format(", ".join(values))

        evaluate_file = self.create_file(os.path.join(ROOT, "generatedtests", ecf["OUTPUT"]), [], "evaluate.py")

        evaluate_file.write("import os\n")
        evaluate_file.write("import sys\n\n")
        evaluate_file.write("sys.path.append(os.path.join(\"{0}\", \"benchmarks\"))\n\n".format(self.ROOT))
        evaluate_file.write("from evaluator import *\n\n")
        evaluate_file.write("if __name__ == '__main__':\n")
        evaluate_file.write("    evaluator = ESPRESOTestEvaluator('{0}')\n".format(os.path.dirname(evaluate_file.name)))

        for table in ecf["TABLES"]:
            rows = ecf["TABLES"][table]["ROWS"].split(",")
            columns = ecf["TABLES"][table]["COLUMNS"].split(",")
            rows = map(lambda x: x.strip(), rows)
            columns = map(lambda x: x.strip(), columns)
            for row in rows:
                if row not in ["VALUES"] + ["L{0}".format(i) for i in range(1, levels + 1)]:
                    print "Table {0} contains invalid ROWS value {1}.".format(table, row)
            for column in columns:
                if column not in ["VALUES"] + ["L{0}".format(i) for i in range(1, levels + 1)]:
                    print "Table {0} contains invalid COLUMNS value {1}.".format(table, column)

            evaluate_file.write("    evaluator.evaluate(\n")
            evaluate_file.write("        tablename='{0}',\n".format(table))
            evaluate_file.write("        rows=['{0}'],\n".format("', '".join(rows)))
            evaluate_file.write("        columns=['{0}'],\n".format("', '".join(columns)))
            for level in range(1, levels + 1):
                if ecf["TABLES"][table]["L{0}".format(level)] == "ALL":
                    values = variables["L{0}".format(level)]
                else:
                    values = self.evaluate(ecf["TABLES"][table]["L{0}".format(level)], variables, "TABLES::{0}::L{1}".format(table, level))
                evaluate_file.write("        L{0}={1},\n".format(level, values))
            values = self.evaluate(ecf["TABLES"][table]["VALUES"], variables, "TABLES::{0}::VALUES".format(table))
            evaluate_file.write("        VALUES={0})\n".format(values))

    def process_file(self, file):
        parser = ECFParser()
        ecf = parser.get_ecf_tree(parser.get_object(parser.tokenize(file), "PYTHON_TEST_GENERATOR"))
        ecf["origin"] = file

        run = self.create_dir(os.path.join(ROOT, "generatedtests", ecf["OUTPUT"], "run"))
        log = self.create_dir(os.path.join(ROOT, "generatedtests", ecf["OUTPUT"], "log"))
        stats = self.create_dir(os.path.join(ROOT, "generatedtests", ecf["OUTPUT"], "stats"))
        self.generate_level(ecf, run, log, stats)
        self.generate_evaluate_file(ecf, log, stats)

    def generate(self, path):
        if os.path.isfile(os.path.join(RUNDIR, path)):
            self.process_file(os.path.join(RUNDIR, path))

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "Provide *.ecf file as the first argument."
        exit()

    generator = ESPRESOTestGenerator();
    generator.generate(sys.argv[1])
