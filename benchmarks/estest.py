
import shutil, os, subprocess, copy
from nose.plugins.skip import SkipTest

try:
    import requests, git
    snailwatch = True
except ImportError:
    snailwatch = False

class ESPRESOTest:

    checker = os.getenv("ECFCHECKER", False)
    fast = os.getenv("FAST", False)
    skip_hypre = os.getenv("SKIP_HYPRE", False)
    skip_mklpdss = os.getenv("SKIP_MKLPDSS", False)
    oversub = os.getenv("OVERSUB", False)
    create = os.getenv("CREATE", False)

    root = os.path.dirname(os.path.dirname(__file__))
    feti4itester = os.path.join(root, "build", "test.feti4i")
    espreso = os.path.join(root, "build", "espreso")
    ecfchecker = os.path.join(root, "build", "ecfchecker")
    waf = os.path.join(root, "waf")

    env = dict(os.environ)

    if "MPIRUN" in env:
        mpirun = env["MPIRUN"].split()
    else:
        mpirun = [ "mpirun", "-n" ]
    if "MPILIB" in env:
        mpilib = env["MPILIB"].strip()
    else:
        mpilib = ""
    if "SYSTEM" in env:
        system = env["SYSTEM"].strip()
    else:
        system = "workstation"
    env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"

    path = ""
    ecf = "espreso.ecf"
    processes = 1
    args = []
    store_results = False
    external = False

    _program = []

    @staticmethod
    def has_snailwatch():
        return snailwatch and "SNAILWATCH_URL" in os.environ and "SNAILWATCH_TOKEN"in os.environ

    def get_commit():
        def is_git_directory(path = '.'):
            return subprocess.call(['git', '-C', path, 'status'], stderr=subprocess.STDOUT, stdout = open(os.devnull, 'w')) == 0

        if is_git_directory():
            return subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip().decode()

    def get_info():
        info = subprocess.check_output([ESPRESOTest.waf, "info"]).rstrip().decode()
        for line in info.splitlines():
            if line.startswith("commit"):
                return line

    @staticmethod
    def parse_run_command(command):
        mpirun = command.split("espreso")[0].split()
        procs = mpirun.index("-n")
        mpirun.pop(procs)
        processes = mpirun[procs]
        mpirun.pop(procs)
        mpirun = mpirun + ["-n"]
        args = command.split("espreso.ecf")[1].split()
        return mpirun, processes, args

    @staticmethod
    def set_threads(threads):
        threads = int(threads)
        ESPRESOTest.env["MKL_NUM_THREADS"] = str(1)
        ESPRESOTest.env["OMP_NUM_THREADS"] = str(threads)
        ESPRESOTest.env["SOLVER_NUM_THREADS"] = str(threads)
        ESPRESOTest.env["PAR_NUM_THREADS"] = str(threads)


    @staticmethod
    def has_solver(solver):
        program = [ ESPRESOTest.waf, "show" ]
        output, error = subprocess.Popen(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        return output.decode().find(solver) != -1

    @staticmethod
    def raise_error(error, output=""):
        raise Exception("\n {4} \n\nPath: {3}\n\nProgram: {2}\n\nERROR:{0}\n\nOUTPUT{1}\n\n {4} \n\n\n".format(
                error, output, " ".join(ESPRESOTest._program), ESPRESOTest.path, "#" * 80))

    @staticmethod
    def clean(path="results"):
        shutil.rmtree(os.path.join(ESPRESOTest.path, path), ignore_errors=True)
        ESPRESOTest.args = []
        ESPRESOTest.store_results = False
        ESPRESOTest.external = False

    @staticmethod
    def run_program(program):
        ESPRESOTest._program = copy.deepcopy(program)
        program.append("--OUTPUT::LOGGER=PARSER")
        program.append("-vvttt")
        if not ESPRESOTest.store_results:
            program.append("--OUTPUT::RESULTS_STORE_FREQUENCY=NEVER")
            program.append("--OUTPUT::MODE=SYNC")
        if ESPRESOTest.has_snailwatch():
            program.append("-m")

        def _popen():
            return subprocess.Popen(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ESPRESOTest.path, env=ESPRESOTest.env)

        p = _popen()
        output, error = p.communicate()

        if p.returncode: # if the program is killed try to run it again
            p = _popen()
            output, error = p.communicate()

        return (output.decode(), error.decode())

    @staticmethod
    def run():
        if "HYPRE" in ESPRESOTest.args and (ESPRESOTest.skip_hypre or not ESPRESOTest.has_solver("hypre")):
            raise SkipTest("Test require installed 'hypre' library.")
        if "MKLPDSS" in ESPRESOTest.args and (ESPRESOTest.skip_mklpdss or not ESPRESOTest.has_solver("mklpdss")):
            raise SkipTest("Test require installed 'mklpdss' library.")
        if ESPRESOTest.external and not os.path.exists(os.path.join(os.sep, "data", "espreso", "mesiotest")):
            raise SkipTest("Test depends on the external file.")
        program = copy.deepcopy(ESPRESOTest.mpirun)
        program.append(str(ESPRESOTest.processes))
        if (ESPRESOTest.oversub):
            program.append("--oversubscribe")
        if ESPRESOTest.checker:
            program.append(ESPRESOTest.ecfchecker)
        else:
            program.append(ESPRESOTest.espreso)
        program.extend([ "-c", os.path.join(ESPRESOTest.path, ESPRESOTest.ecf) ])
        program.extend(map(str, ESPRESOTest.args))

        output, error = ESPRESOTest.run_program(program)
        if error != "":
            ESPRESOTest.raise_error(error, output)
        return output

    @staticmethod
    def compare_emr(preset):
        if ESPRESOTest.checker:
            return

        def compare(v1, v2):
            if v1 != v2 and abs(float(v1) - float(v2)) > 1e-4 and abs((float(v1) - float(v2)) / float(v1)) > 1e-3:
                error = "various monitored results:\n"
                error += "  region={0}".format(table1[0][column])
                error += ", property={0}".format(table1[1][column])
                error += ", stat={0}".format(table1[2][column])
                error += ", step={0}".format(table1[row][0])
                error += ", substep={0}\n".format(table1[row][1])
                error += "  preset={0} != current={1}\n".format(v1, v2)
                ESPRESOTest.raise_error(error)

        def create_table(emr):
            return [ [ value.strip() for value in line.split(";") ] for line in open(emr, "r").readlines() if len(line.strip()) ]

        emr = os.path.join(ESPRESOTest.path, "results", "last", ESPRESOTest.ecf.replace(".ecf", ".emr"))
        if not os.path.isfile(emr):
            ESPRESOTest.raise_error("Missing monitoring report '{0}'.".format(emr))

        table1 = create_table(os.path.join(ESPRESOTest.path, preset))
        table2 = create_table(emr)

        if len(table1) != len(table2):
            ESPRESOTest.raise_error("various time steps")
        if len(table1[0]) != len(table2[1]):
            ESPRESOTest.raise_error("various monitored properties")

        for row, (row1, row2) in enumerate(zip(table1, table2)):
            for column, (value1, value2) in enumerate(zip(row1, row2)):
                compare(value1.strip(), value2.strip())

    @staticmethod
    def create_emr(preset):
        if ESPRESOTest.checker:
            return
        
        emr = os.path.join(ESPRESOTest.path, "results", "last", ESPRESOTest.ecf.replace(".ecf", ".emr"))
        if not os.path.isfile(emr):
            ESPRESOTest.raise_error("Missing monitoring report '{0}'.".format(emr))

        try:
            shutil.copyfile(emr, os.path.join(ESPRESOTest.path, preset))
        except OSError:
            ESPRESOTest.raise_error("Cannot copy emr file")

    @staticmethod
    def compare_mesh(preset, output):
        if ESPRESOTest.checker:
            return

        def read_log(log):
            mesh = dict();
            lines = [ line.split(": ")[1].strip() for line in log if line.startswith("mesh:") ]
            for line in lines:
                for param in line.split(", "):
                    key, value = param.split("=")
                    if key == "region":
                        region = value
                        if region not in mesh:
                            mesh[value] = dict()
                    else:
                        mesh[region][key] = value
            return mesh

        log1 = read_log(open(os.path.join(ESPRESOTest.path, preset), "r").readlines())
        log2 = read_log(output.splitlines())

        for region, data in log1.items():
            for key, value in data.items():
                if region not in log2:
                    ESPRESOTest.raise_error("missing mesh region: {0}\n".format(region), output)
                if key not in log2[region]:
                    ESPRESOTest.raise_error("missing key='{0}' on mesh region: {1}\n".format(key, region), output)
                if log2[region][key] != value:
                    ESPRESOTest.raise_error("invalid value on region {0}::{1} {2} != {3}\n".format(region, key, value, log2[region][key]), output)

    @staticmethod
    def check_continuity(output):
        if ESPRESOTest.checker:
            return

        for line in output.splitlines():
            if line.startswith("decomposition: "):
                if int(line.split(":")[2]) != 1:
                    ESPRESOTest.raise_error("non-continuous cluster found.\n")

    @staticmethod
    def extract(output, values):
        if ESPRESOTest.checker:
            return

        def parse(line):
            stats = dict()
            for stat in line.split(",")[1:]:
                stats[str(stat.split("=")[0].strip(" '"))] = float(stat.split("=")[1].strip(" '"))
            return stats

        measurement = dict()
        for line in output.splitlines():
            for value in values:
                if value in line:
                    measurement[value] = parse(line)
                    break
        return measurement

    @staticmethod
    def report_mesurement(tab, functions):
        output = open(os.path.join(ESPRESOTest.path, ESPRESOTest.get_info().lstrip("commit-") + ".csv"), "w")
        output.write("{0:{width}}".format("", width=50))
        for instance in tab:
            output.write(",{0:>{width}}".format("_".join(map(str, instance[0])), width=15))
        output.write("\n")
        for fnc in functions:
            output.write("{0:{width}}".format("{} [{}]".format(fnc[0], fnc[1].upper()), width=50))
            for instance in tab:
                output.write(",{0:>{width}}".format(instance[1][fnc[0]][fnc[1]], width=15))
            output.write("\n")

    @staticmethod
    def report(timereport):
        if not ESPRESOTest.has_snailwatch() or ESPRESOTest.checker:
            return

        header = {
            "Content-Type": "application/json",
            "Authorization": os.environ['SNAILWATCH_TOKEN']
        }

        path = os.path.relpath(ESPRESOTest.path, os.path.dirname(__file__))

        benchmark = "-".join(path.split("/") + map(str, ESPRESOTest.args))
        commit = git.Repo(search_parent_directories=True).head.object.hexsha
        processes = str(ESPRESOTest.processes)
        threads = ESPRESOTest.env["OMP_NUM_THREADS"]

        log = os.path.join(ESPRESOTest.path, "results", "last", ESPRESOTest.ecf.replace(".ecf", ".log"))

        results = { }
        regions = ["Mesh preprocessing timing- Total", "Physics solver timing- Total"]
        with open(log, 'r') as file:
            for line in file:
                for region in regions:
                    if line.startswith(region):
                        results[region.replace(" ", "_")] = { "type": "time", "value": line.split("avg.:")[1].split()[0] }

        response = requests.post(
            os.environ['SNAILWATCH_URL'],
            headers=header,
            json={
                "benchmark": benchmark,
                "environment": { "commit": commit, "processes": processes, "threads": threads },
                "result": results
            })

        if response.status_code != 201:
            ESPRESOTest.raise_error("Cannot push to snailwatch")

