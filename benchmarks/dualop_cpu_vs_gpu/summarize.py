#!/usr/bin/env python3

import os
import sys
import csv
from datetime import datetime
import io

def read_csv_to_arrays(file_path):
    rows = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file, delimiter=";")
        rows = [row for row in csv_reader]
        return rows

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







outfile = summ_dir + "/timings.csv"
outstring = io.StringIO()
cells_info = ["id", "str", "run_id", "machine", "tool", "dual_operator", "problem", "physics", "dim", "dofs_per_node", "element_type", "domains_x", "domains_y", "domains_z", "n_domains_total", "elements_x", "elements_y", "elements_z", "n_elements_per_domain", "avg_domain_surface", "avg_domain_volume", "n_dofs"]
cells_options = ["CONCURRENCY_SET", "CONCURRENCY_UPDATE", "TRS1_FACTOR_STORAGE", "TRS1_SOLVE_TYPE", "TRS2_FACTOR_STORAGE", "TRS2_SOLVE_TYPE", "TRSM_RHS_SOL_ORDER", "PATH_IF_HERMITIAN", "CONCURRENCY_APPLY", "APPLY_SCATTER_GATHER_WHERE"]
outstring.write(";".join(cells_info))
outstring.write(";;")
outstring.write("error;timeout")
outstring.write(";;")
outstring.write("set;update;apply;;set/subdomain;update/subdomain;apply/subdomain")
outstring.write(";;")
outstring.write(";".join(cells_options))
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
        set_time_perdom = str(float(set_time) / n_domains_total) if len(set_time)>0 else ""
        update_time_perdom = str(float(update_time) / n_domains_total) if len(update_time)>0 else ""
        apply_avg_perdom = apply_avg / n_domains_total
        apply_avg_perdom_str = str(apply_avg_perdom) if len(apply_lines) > 0 else ""

        outstring.write(set_time + ";")
        outstring.write(update_time + ";")
        outstring.write(apply_avg_str + ";")
        outstring.write(";")
        outstring.write(set_time_perdom + ";")
        outstring.write(update_time_perdom + ";")
        outstring.write(apply_avg_perdom_str + ";")
        outstring.write(";")

        for opt in cells_options:
            lines = list(filter(lambda line: opt in line, output_lines))
            if len(lines) > 0:
                outstring.write(lines[0][70:-2].replace(" ", ""))
            outstring.write(";")

        outstring.write("\n")
write_string_to_file(outfile, outstring.getvalue())
outstring.close()






















machines_tools_dualops = [
    ("karolina", "cudalegacy",   "EXPLICIT_GPU"),
    ("karolina", "cudalegacy",   "IMPLICIT_GPU"),
    ("karolina", "cudamodern",   "EXPLICIT_GPU"),
    ("karolina", "cudamodern",   "IMPLICIT_GPU"),
    ("karolina", "mklpardiso",   "EXPLICIT_SC"),
    ("karolina", "mklpardiso",   "IMPLICIT"),
    ("karolina", "suitesparse",  "EXPLICIT_SC"),
    ("karolina", "suitesparse",  "IMPLICIT"),
    ("lumi",     "rocm",         "EXPLICIT_GPU"),
    ("lumi",     "rocm",         "IMPLICIT_GPU"),
    ("lumi",     "mklpardiso",   "EXPLICIT_SC"),
    ("lumi",     "mklpardiso",   "IMPLICIT"),
    ("lumi",     "suitesparse",  "EXPLICIT_SC"),
    ("lumi",     "suitesparse",  "IMPLICIT"),
    ("e4red",    "cudamodern",   "EXPLICIT_GPU"),
    ("e4red",    "cudamodern",   "IMPLICIT_GPU"),
    ("e4red",    "suitesparse",  "EXPLICIT_SC"),
    ("e4red",    "suitesparse",  "IMPLICIT"),
    ("tiber",    "oneapi",       "EXPLICIT_GPU"),
    ("tiber",    "oneapi",       "IMPLICIT_GPU"),
    ("tiber",    "mklpardiso",   "EXPLICIT_SC"),
    ("tiber",    "mklpardiso",   "IMPLICIT"),
    ("tiber",    "suitesparse",  "EXPLICIT_SC"),
    ("tiber",    "suitesparse",  "IMPLICIT"),
    ("sprddr",   "mklpardiso",   "EXPLICIT_SC"),
    ("sprddr",   "mklpardiso",   "IMPLICIT"),
    ("sprddr",   "suitesparse",  "EXPLICIT_SC"),
    ("sprddr",   "suitesparse",  "IMPLICIT"),
    ("sprhbm",   "mklpardiso",   "EXPLICIT_SC"),
    ("sprhbm",   "mklpardiso",   "IMPLICIT"),
    ("sprhbm",   "suitesparse",  "EXPLICIT_SC"),
    ("sprhbm",   "suitesparse",  "IMPLICIT")
]

csv_all = read_csv_to_arrays(summ_dir + "/timings.csv")
csv_header = csv_all[0]
csv_data = csv_all[1:]

machine_col = csv_header.index("machine")
tool_col = csv_header.index("tool")
dualop_col = csv_header.index("dual_operator")
problem_col = csv_header.index("problem")
element_type_col = csv_header.index("element_type")
n_domains_col = csv_header.index("n_domains_total")
n_dofs_col = csv_header.index("n_dofs")
set_col = csv_header.index("set")
update_col = csv_header.index("update")
apply_col = csv_header.index("apply")
set_ps_col = csv_header.index("set/subdomain")
update_ps_col = csv_header.index("update/subdomain")
apply_ps_col = csv_header.index("apply/subdomain")

outfile = summ_dir + "/compare.csv"
outstring = io.StringIO()
cells_info = ["id", "str", "run_id", "problem", "physics", "dim", "dofs_per_node", "element_type", "elements_x", "elements_y", "elements_z", "n_elements_per_domain", "avg_domain_surface", "avg_domain_volume", "n_dofs"]
cells_info_cols = [csv_header.index(ci) for ci in cells_info]
cells_mtd = ["machine1", "tool1", "dual_operator1", "machine2", "tool2", "dual_operator2"]
outstring.write(";".join(cells_info))
outstring.write(";;")
outstring.write(";".join(cells_mtd))
outstring.write(";;")
outstring.write("ndom1;ndom2;;set1;set2;update1;update2;apply1;apply2;;set1/subdomain;set2/subdomain;update1/subdomain;update2/subdomain;apply1/subdomain;apply2/subdomain")
outstring.write("\n")



n_td = len(machines_tools_dualops)
for i in range(n_td):
    machine_i = machines_tools_dualops[i][0]
    tool_i = machines_tools_dualops[i][1]
    dualop_i = machines_tools_dualops[i][2]
    csv_data_i = list(filter(lambda row: row[machine_col] == machine_i and row[tool_col] == tool_i and row[dualop_col] == dualop_i, csv_data))

    for j in range(n_td):
        machine_j = machines_tools_dualops[j][0]
        tool_j = machines_tools_dualops[j][1]
        dualop_j = machines_tools_dualops[j][2]
        csv_data_j = list(filter(lambda row: row[machine_col] == machine_j and row[tool_col] == tool_j and row[dualop_col] == dualop_j, csv_data))

        for row_i in csv_data_i:
            rows_j = list(filter(lambda row: row[problem_col] == row_i[problem_col] and row[element_type_col] == row_i[element_type_col] and row[n_dofs_col] == row_i[n_dofs_col], csv_data_j))
            if len(rows_j) != 1: print("ERROR too many rows")
            row_j = rows_j[0]

            for ci_idx in range(len(cells_info)):
                from_i = row_i[cells_info_cols[ci_idx]]
                from_j = row_j[cells_info_cols[ci_idx]]
                outstring.write((from_i if len(from_i)>0 else from_j).replace("TETRA4", "TETRA04") + ";")
            outstring.write(";")
            outstring.write(machine_i + ";")
            outstring.write(tool_i + ";")
            outstring.write(dualop_i + ";")
            outstring.write(machine_j + ";")
            outstring.write(tool_j + ";")
            outstring.write(dualop_j + ";")
            outstring.write(";")
            outstring.write(row_i[n_domains_col] + ";")
            outstring.write(row_j[n_domains_col] + ";")
            outstring.write(";")
            outstring.write(row_i[set_col] + ";")
            outstring.write(row_j[set_col] + ";")
            outstring.write(row_i[update_col] + ";")
            outstring.write(row_j[update_col] + ";")
            outstring.write(row_i[apply_col] + ";")
            outstring.write(row_j[apply_col] + ";")
            outstring.write(";")
            outstring.write(row_i[set_ps_col] + ";")
            outstring.write(row_j[set_ps_col] + ";")
            outstring.write(row_i[update_ps_col] + ";")
            outstring.write(row_j[update_ps_col] + ";")
            outstring.write(row_i[apply_ps_col] + ";")
            outstring.write(row_j[apply_ps_col] + ";")
            outstring.write("\n")
write_string_to_file(outfile, outstring.getvalue())
outstring.close()
