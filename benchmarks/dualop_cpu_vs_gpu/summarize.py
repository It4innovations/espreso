#!/usr/bin/env python3

import os
import sys
from datetime import datetime
import io

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

def get_n_dofs_per_node(problem):
    if problem == "heat_transfer_2d":
        return 1
    elif problem == "heat_transfer_3d":
        return 1
    elif problem == "linear_elasticity_2d":
        return 2
    elif problem == "linear_elasticity_3d":
        return 3
    else:
        return -1

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




basedir = "benchmarks/dualop_cpu_vs_gpu"
if not os.path.exists(basedir + "/summarize.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)



# copy or put symlinks to this folder
results_dir = basedir + "/results_to_summarize"

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")
summ_dir = basedir + "/summary/" + datestr
os.makedirs(summ_dir, exist_ok=True)







outfile = summ_dir + "/compare.csv"
outstring = io.StringIO()
cells_info = ["id", "str", "run_id", "machine", "tool", "dual_operator", "problem", "physics", "dim", "dofs_per_node", "element_type", "domains_x", "domains_y", "domains_z", "n_domains_total", "elements_x", "elements_y", "elements_z", "n_elements_per_domain", "avg_domain_surface", "avg_domain_volume", "n_dofs"]
outstring.write(";".join(cells_info))
outstring.write(";;")
outstring.write("error;timeout")
outstring.write(";;")
outstring.write("set;update;apply")
outstring.write("\n")
for dir_name in os.listdir(results_dir):
    dir_path = results_dir + "/" + dir_name
    for run_id in os.listdir(dir_path):
        run_path = dir_path + "/" + run_id
        info_lines = read_file_to_string(run_path + "/info.txt").split("\n")
        problem = list(filter(lambda line: "ecf_file" in line, info_lines))[0].split(".")[1]
        output_lines = read_file_to_string(run_path + "/last/espreso." + problem + ".log").split("\n")
        timer_lines = list(filter(lambda line: "TMP DUAL OPERATOR" in line, output_lines))
        n_domains_total = 1
        n_elements_per_domain = 1
        n_dofs_per_domain = 1
        elem_type = ""
        for field in cells_info:
            if field == "id" or field == "str":
                val = ""
            elif field == "physics":
                val = problem[:-3]
            elif field == "dim":
                val = problem[-2:-1]
            elif field == "n_domains_total":
                val = str(n_domains_total)
            elif field == "n_elements_per_domain":
                val = str(n_elements_per_domain)
            elif field == "n_dofs":
                val = str(n_dofs_per_domain * get_n_dofs_per_node(problem))
            elif field == "run_id":
                val = run_id
            elif field == "problem":
                val = problem
            elif field == "dofs_per_node":
                val = str(get_n_dofs_per_node(problem))
            elif field == "avg_domain_volume":
                lines = list(filter(lambda line: "Domain volume [dofs]" in line, output_lines))
                if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                else:              val = ""
            elif field == "avg_domain_surface":
                lines = list(filter(lambda line: "Domain surface [dofs]" in line, output_lines))
                if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                else:              val = ""
            else:
                val = list(filter(lambda line: field in line, info_lines))[0].split(" ")[1]
                if field == "element_type":
                    elem_type = val
                if field.startswith("domains_"):
                    n_domains_total = n_domains_total * int(val)
                if field.startswith("elements_"):
                    n_elements_per_domain = n_elements_per_domain * int(val)
                    if not (problem[-2:-1] == "2" and field == "elements_z"):
                        n_dofs_per_domain = n_dofs_per_domain * get_nnodes_per_side(int(val), elem_type)
            outstring.write(val)
            outstring.write(";")

        outstring.write(";")
        outstring.write("\"" + read_file_to_string(run_path + "/last/espreso." + problem + ".err").split("\0")[0].replace("\n", ";").replace("\"", "\"\"") + "\";")
        outstring.write("\"" + read_file_to_string(run_path + "/timeout.txt").replace("\n", ";").replace("\"", "\"\"") + "\";")
        outstring.write(";")

        set_lines = list(filter(lambda line: "TMP DUAL OPERATOR SET TIME:" in line, timer_lines))
        set_time = set_lines[0][31:43].replace(" ", "") if len(set_lines) > 0 else ""
        update_lines = list(filter(lambda line: "TMP DUAL OPERATOR UPDATE TIME:" in line, timer_lines))
        update_time = update_lines[-1][31:43].replace(" ", "") if len(update_lines) > 0 else ""
        apply_lines = list(filter(lambda line: "TMP DUAL OPERATOR APPLY TIME:" in line, timer_lines))
        apply_times = [float(line[31:43]) for line in apply_lines[-10:]] if len(apply_lines) > 0 else []
        apply_avg = sum(apply_times) / len(apply_times) if len(apply_lines) > 0 else 0
        apply_avg_str = str(apply_avg) if len(apply_lines) > 0 else ""
        outstring.write(set_time + ";")
        outstring.write(update_time + ";")
        outstring.write(apply_avg_str)

        outstring.write("\n")
write_string_to_file(outfile, outstring.getvalue())
outstring.close()
