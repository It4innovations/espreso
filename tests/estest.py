
import shutil, os, subprocess, copy, re

class ESPRESOTest:

    root = os.path.dirname(os.path.dirname(__file__))
    espreso = os.path.join(root, "build", "espreso")
    env = dict(os.environ)

    mpirun = [ "mpirun", "-n" ]
    env["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    env["SOLVER_NUM_THREADS"] = "1"
    env["PAR_NUM_THREADS"] = "1"
    env["OMPI_MCA_rmaps_base_mapping_policy"] = "core"

    path = ""
    ecf = "espreso.ecf"
    processes = 1
    args = []
    store_results = False
    external = False

    _program = []

    def get_info():
        info = subprocess.check_output([ESPRESOTest.waf, "info"]).rstrip().decode()
        for line in info.splitlines():
            if line.startswith("commit"):
                return line

    @staticmethod
    def set_threads(threads):
        threads = int(threads)
        ESPRESOTest.env["MKL_NUM_THREADS"] = str(1)
        ESPRESOTest.env["OMP_NUM_THREADS"] = str(threads)
        ESPRESOTest.env["SOLVER_NUM_THREADS"] = str(threads)
        ESPRESOTest.env["PAR_NUM_THREADS"] = str(threads)
        ESPRESOTest.env["OMPI_MCA_rmaps_base_mapping_policy"] = "slot:pe=" + str(threads)

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
        # program.append("--LOOP=CONDITIONS")
        program.append("-ttt")
        if not ESPRESOTest.store_results:
            program.append("--OUTPUT::RESULTS_STORE_FREQUENCY=NEVER")
            program.append("--OUTPUT::MODE=SYNC")

        p = subprocess.Popen(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ESPRESOTest.path, env=ESPRESOTest.env)
        output, error = p.communicate()
        return (output.decode(), error.decode())

    @staticmethod
    def run():
        program = copy.deepcopy(ESPRESOTest.mpirun)
        program.append(str(ESPRESOTest.processes))
        program.append(ESPRESOTest.espreso)
        program.extend([ "-c", os.path.join(ESPRESOTest.path, ESPRESOTest.ecf) ])
        program.extend(map(str, ESPRESOTest.args))

        output, error = ESPRESOTest.run_program(program)
        if error != "":
            ESPRESOTest.raise_error(error, output)
        return output

    @staticmethod
    def compare_emr(preset):
        def compare(v1, v2):
            if v1 != v2 and abs(float(v1) - float(v2)) > 1e-4 and abs((float(v1) - float(v2)) / float(v2)) > 1e-3:
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

        preset = preset.strip()
        emr = os.path.join(ESPRESOTest.path, "results", "last", ESPRESOTest.ecf.replace(".ecf", ".emr"))
        if not os.path.isfile(emr):
            ESPRESOTest.raise_error("Missing monitoring report '{0}'.".format(emr))

        if not os.path.isfile(os.path.join(ESPRESOTest.path, preset)):
            shutil.copy(emr, os.path.join(ESPRESOTest.path, preset))
            ESPRESOTest.raise_error("Monitoring report '{0}' created.".format(preset))

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
    def compare_mesh(preset, output):
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
    def extract(output, values):
        def parse(key, line):
            stats = dict()
            if line.startswith("time"):
                for stat in line.split(",")[1:]:
                    stats[str(stat.split("=")[0].strip(" '"))] = float(stat.split("=")[1].strip(" '"))
                return stats
            else:
                return line.replace(key, "").strip(" =");

        measurement = dict()
        for line in output.splitlines():
            for value in values:
                if value in line:
                    measurement[value] = parse(value, line)
                    break
        return measurement

    @staticmethod
    def parse(tab, key, stats):
        for header, line in tab:
            if key in line and isinstance(line[key], str):
                values = []
                for v in line[key].split():
                    try:
                        v = v.strip(" ><")
                        if "." in v:
                            value = float(v)
                        else:
                            value = int(v)
                        values.append(value)
                    except ValueError:
                        pass
                line[key] = dict()
                for i, value in enumerate(values):
                    line[key][stats[i]] = value

    @staticmethod
    def report_mesurement(name, tab, functions):
        output = open(os.path.join(ESPRESOTest.path, name + "." + ESPRESOTest.get_info().lstrip("commit-") + ".csv"), "w")
        output.write("{0:{width}}".format("", width=50))
        size = 0
        for instance in tab:
            size = max(size, len("_".join(map(str, instance[0]))) + 1)
        for instance in tab:
            output.write(",{0:>{width}}".format("_".join(map(str, instance[0])), width=size))
        output.write("\n")
        for fnc in functions:
            output.write("{0:{width}}".format("{} [{}]".format(fnc[0], fnc[1].upper()), width=50))
            for instance in tab:
                if fnc[0] in instance[1]:
                    output.write(",{0:>{width}}".format(instance[1][fnc[0]][fnc[1]], width=size))
                else:
                    output.write(",{0:>{width}}".format("", width=size))
            output.write("\n")

