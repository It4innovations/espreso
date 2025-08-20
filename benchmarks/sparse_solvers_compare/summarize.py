#!/usr/bin/env python3

import os
import sys
import io
from datetime import datetime



def read_file_to_string(file_path):
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        print("could not open file " + file_path)
        sys.exit(2)

def write_string_to_file(file_path, str):
    try:
        with open(file_path, 'w') as file:
            file.write(str)
    except FileNotFoundError:
        print("could not open file " + file_path)
        sys.exit(3)

def get_nnodes_per_side(nelems, element_type):
    if element_type == "TRIANGLE3":
        return nelems + 1
    if element_type == "TRIANGLE6":
        return 2 * nelems + 1
    if element_type == "TETRA4":
        return nelems + 1
    if element_type == "TETRA10":
        return 2 * nelems + 1
    return 0



if len(sys.argv) <= 1:
    print("not enough arguments")
    exit(1)
run_dir_name = sys.argv[1]

basedir = "benchmarks/sparse_solvers_compare"
run_dir = basedir + "/runs/" + run_dir_name

summdir = run_dir + "/summary"
os.makedirs(summdir, exist_ok=True)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")
summfile = summdir + "/" + datestr + ".csv"

summstring = io.StringIO()





header_1 = ["my_id", "task_id", "physics", "dimension"]
header_2 = ["element_type", "domains_x", "domains_y", "domains_z", "domains_total", "elements_x", "elements_y", "elements_z", "elements_total", "n_dofs", "n_dofs_domain_actual", "n_dofs_interface_actual"]
header_3 = ["dualoperator", "schur_impl", "sparse_solver_impl", "explicit_apply_where"]
header_4 = ["exitcode", "errfile", "timeoutfile"]
header_5 = ["update_time", "apply_time", "update_time_per_subdomain", "apply_time_per_subdomain"]
summstring.write(";".join(header_1))
summstring.write(";;")
summstring.write(";".join(header_2))
summstring.write(";;")
summstring.write(";".join(header_3))
summstring.write(";;")
summstring.write(";".join(header_4))
summstring.write(";;")
summstring.write(";".join(header_5))
summstring.write("\n")



tasks_dir = run_dir + "/tasks"
for task_dir_name in os.listdir(tasks_dir):
    task_dir = tasks_dir + "/" + task_dir_name
    task_id = task_dir_name
    print(task_dir)

    config_file_path = task_dir + "/config.sh"
    config_file_lines = read_file_to_string(config_file_path).split("\n")

    ecf_line = [line for line in config_file_lines if ("ESPRESO_BENCHMARK_ecf_file" in line)][0]
    problem = ecf_line.split(".")[1]
    physics = problem[0:-3]
    dimension = problem[-2:-1]

    summstring.write(";")
    summstring.write(task_id + ";")
    summstring.write(physics + ";")
    summstring.write(dimension + ";")
    summstring.write(";")
    
    out_file = task_dir + "/log_out.txt"
    out_file_lines = read_file_to_string(out_file).split("\n")

    domains_total = 1
    elements_total = 1
    n_dofs = 1
    element_type = ""
    for field in header_2:
        if field == "":
            summstring.write(";")
        elif field == "domains_total":
            summstring.write(str(domains_total) + ";")
        elif field == "elements_total":
            summstring.write(str(elements_total) + ";")
        elif field == "n_dofs":
            summstring.write(str(n_dofs) + ";")
        elif field == "n_dofs_domain_actual":
            lines = [l for l in out_file_lines if (" =   Domain volume [dofs]" in l)]
            if len(lines) == 0:
                summstring.write(";")
            else:
                line = lines[0]
                ndofs = line.split("[dofs]")[1].split("<")[0].replace(" ", "")
                summstring.write(ndofs + ";")
        elif field == "n_dofs_interface_actual":
            lines = [l for l in out_file_lines if (" =   Domain surface [dofs]" in l)]
            if len(lines) == 0:
                summstring.write(";")
            else:
                line = lines[0]
                ndofs = line.split("[dofs]")[1].split("<")[0].replace(" ", "")
                summstring.write(ndofs + ";")
        else:
            line = [line for line in config_file_lines if (("ECF_" + field.upper()) in line)][0]
            value = line.split("=")[1]
            summstring.write(value + ";")
            if field == "element_type":
                element_type = value
            if "domains" in field:
                domains_total *= int(value)
            if "elements" in field:
                elements_total *= int(value)
                if not (dimension == "2" and field == "elements_z"):
                    n_dofs *= get_nnodes_per_side(int(value), element_type)
    summstring.write(";")
    for field in header_3:
        if field == "":
            summstring.write(";")
        else:
            line = [line for line in config_file_lines if (("ECF_" + field.upper()) in line)][0]
            value = line.split("=")[1]
            summstring.write(value + ";")
    summstring.write(";")

    finished_file = task_dir + "/.finished"
    exitcode = read_file_to_string(finished_file).replace("\n", "")
    summstring.write(exitcode + ";")

    err_file = task_dir + "/stderr.txt"
    err_content = read_file_to_string(err_file).replace("\n", ";").replace("\"", "\"\"")
    summstring.write("\"" + err_content + "\"" + ";")
    
    timeout_file = task_dir + "/timeout.txt"
    timeout_content = read_file_to_string(timeout_file).replace("\n", ";").replace("\"", "\"\"")
    summstring.write("\"" + timeout_content + "\"" + ";")

    summstring.write(";")

    measured_lines_update = [line for line in out_file_lines if ("finishd 'update_mainloop'" in line or "finishd 'update_factorize_numeric'" in line)]
    measured_lines_apply = [line for line in out_file_lines if ("finishd 'dualop_explicit_applicator::apply'" in line or "finishd 'cpu_implicit_apply_compute'" in line)]

    update_time = 0
    if len(measured_lines_update) > 0:
        measured_line = measured_lines_update[-1]
        update_time = float(measured_line[100:-3].replace(" ", ""))
        summstring.write(str(update_time) + ";")
    else:
        summstring.write(";")
    
    apply_time = 0
    if len(measured_lines_apply) > 2:
        measured_lines = measured_lines_apply[2:]
        apply_times = [float(l[100:-3].replace(" ", "")) for l in measured_lines]
        apply_time = sum(apply_times) / len(apply_times)
        summstring.write(str(apply_time) + ";")
    else:
        summstring.write(";")

    if len(measured_lines_update) > 0:
        avg_update_time = update_time / domains_total
        summstring.write(str(avg_update_time) + ";")
    else:
        summstring.write(";")

    if len(measured_lines_apply) > 0:
        avg_apply_time = apply_time / domains_total
        summstring.write(str(avg_apply_time) + ";")
    else:
        summstring.write(";")

    summstring.write("\n")

write_string_to_file(summfile, summstring.getvalue())
summstring.close()



print("finished. later graph with")
print("    ./benchmarks/sparse_solvers_compare/graphs.py " + summfile)
