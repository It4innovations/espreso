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

def get_n_nodes(elements_x, elements_y, elements_z, element_type):
    if element_type == "TRIANGLE3":
        return (elements_x + 1) * (elements_y + 1)
    if element_type == "TRIANGLE6":
        return (2 * elements_x + 1) * (2 * elements_y + 1)
    if element_type == "SQUARE4":
        return (elements_x + 1) * (elements_y + 1)
    if element_type == "SQUARE8":
        return (2 * elements_x + 1) * (2 * elements_y + 1) - (elements_x * elements_y)
    if element_type == "TETRA4":
        return (elements_x + 1) * (elements_y + 1) * (elements_z + 1)
    if element_type == "TETRA10":
        return (2 * elements_x + 1) * (2 * elements_y + 1) * (2 * elements_z + 1)
    if element_type == "HEXA8":
        return (elements_x + 1) * (elements_y + 1) * (elements_z + 1)
    if element_type == "HEXA20":
        return (2 * elements_x + 1) * (2 * elements_y + 1) * (2 * elements_z + 1) - (elements_x * elements_y * elements_z)
    raise ValueError("wrong element type")

def get_n_dofs_per_node(physics, dimension):
    if physics == "heat_transfer":
        return 1
    if physics == "linear_elasticity":
        if dimension == "2":
            return 2
        if dimension == "3":
            return 3
    raise ValueError("wrong arguments")



if len(sys.argv) <= 1:
    print("not enough arguments")
    exit(1)
run_id = sys.argv[1]

basedir = "benchmarks/dualop_cluster_domain_size"
run_dir = basedir + "/runs/" + run_id

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")

summdir = basedir + "/summary/" + datestr
os.makedirs(summdir, exist_ok=True)

summfile = summdir + "/summary.csv"

summstring = io.StringIO()





header_1 = ["my_id", "task_id", "physics", "dimension"]
header_2 = ["feti_method", "dual_operator", "preconditioner", "projector", "element_type", "domains_x", "domains_y", "domains_z", "domains_per_cluster", "elements_x", "elements_y", "elements_z", "elements_per_domain", "n_dofs_domain_calc", "n_dofs_domain_read", "n_dofs_interface_read", "n_dofs_total_calc", "n_dofs_total_read"]
header_4 = ["exitcode", "errfile", "timeoutfile"]
header_5 = ["update_time", "apply_time"]
summstring.write(";".join(header_1))
summstring.write(";;")
summstring.write(";".join(header_2))
summstring.write(";;")
summstring.write(";".join(header_4))
summstring.write(";;")
summstring.write(";".join(header_5))
summstring.write("\n")

index = 0

tasks_dir = run_dir + "/tasks"
for task_dir_name in os.listdir(tasks_dir):
    index = index + 1
    if index % 100 == 0:
        print(index)
    task_dir = tasks_dir + "/" + task_dir_name
    task_id = task_dir_name

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

    domains_x = 1
    domains_y = 1
    domains_z = 1
    elements_x = 1
    elements_y = 1
    elements_z = 1
    n_dofs = 1
    element_type = ""
    n_dofs_per_domain_calc = 0
    for field in header_2:
        if field == "":
            summstring.write(";")
        elif field == "domains_per_cluster":
            summstring.write(str(domains_x * domains_y * domains_z) + ";")
        elif field == "elements_per_domain":
            summstring.write(str(elements_x * elements_y * elements_z) + ";")
        elif field == "n_dofs_domain_calc":
            n_nodes_domain = get_n_nodes(elements_x, elements_y, elements_z, element_type)
            n_dofs_per_node = get_n_dofs_per_node(physics, dimension)
            n_dofs_per_domain_calc = n_nodes_domain * n_dofs_per_node
            summstring.write(str(n_dofs_per_domain_calc) + ";")
        elif field == "n_dofs_domain_read":
            lines = [l for l in out_file_lines if (" =   Domain volume [dofs]" in l)]
            if len(lines) == 0:
                summstring.write(";")
            else:
                line = lines[0]
                ndofs = line.split("[dofs]")[1].split("<")[0].replace(" ", "")
                summstring.write(ndofs + ";")
        elif field == "n_dofs_interface_read":
            lines = [l for l in out_file_lines if (" =   Domain interface [dofs]" in l)]
            if len(lines) == 0:
                summstring.write(";")
            else:
                line = lines[0]
                ndofs = line.split("[dofs]")[1].split("<")[0].replace(" ", "")
                summstring.write(ndofs + ";")
        elif field == "n_dofs_total_calc":
            total_elements_x = elements_x * domains_x
            total_elements_y = elements_y * domains_y
            total_elements_z = elements_z * domains_z
            total_nodes = get_n_nodes(total_elements_x, total_elements_y, total_elements_z, element_type)
            total_dofs = total_nodes * get_n_dofs_per_node(physics, dimension)
            summstring.write(str(total_dofs) + ";")
        elif field == "n_dofs_total_read":
            lines = [l for l in out_file_lines if ("  NODES                    :" in l)]
            if len(lines) == 0:
                summstring.write(";")
            else:
                val_str = lines[0][55:67].replace(",", "").replace(" ", "")
                val = int(val_str) * get_n_dofs_per_node(physics, dimension)
                summstring.write(str(val) + ";")
        else:
            line = [line for line in config_file_lines if (("ESP_" + field.upper()) in line)][0]
            value = line.split("=")[1]
            summstring.write(value + ";")
            if field == "element_type":
                element_type = value
            if "domains" in field:
                if field == "domains_x": domains_x = int(value)
                if field == "domains_y": domains_y = int(value)
                if field == "domains_z": domains_z = int(value)
            if "elements" in field:
                if field == "elements_x": elements_x = int(value)
                if field == "elements_y": elements_y = int(value)
                if field == "elements_z": elements_z = int(value)
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

    measured_lines_update = [l for l in out_file_lines if (("::update'" in l) and ("finishd '" in l))]
    measured_lines_apply = [l for l in out_file_lines if (("::apply (vector)'" in l) and ("finishd '" in l))]

    if len(measured_lines_update) > 0:
        measured_line = measured_lines_update[-1]
        measured_time = measured_line[100:-3].replace(" ", "")
        summstring.write(measured_time + ";")
    else:
        summstring.write(";")

    if len(measured_lines_apply) > 0:
        if len(measured_lines_apply) == 1:
            measured_lines_apply_to_use = measured_lines_apply
        elif len(measured_lines_apply) <= 3:
            measured_lines_apply_to_use = measured_lines_apply[1:]
        else:
            measured_lines_apply_to_use = measured_lines_apply[2:]
        vals = [float(l[100:-3].replace(" ", "")) for l in measured_lines_apply_to_use]
        avg = sum(vals) / len(vals)
        summstring.write(str(avg) + ";")
    else:
        summstring.write(";")

    summstring.write("\n")

write_string_to_file(summfile, summstring.getvalue())
summstring.close()

print("later graph with:")
print("  ./benchmarks/dualop_cluster_domain_size/graphs.py " + summfile)
