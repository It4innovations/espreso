#!/usr/bin/env python3

import os
import sys
import csv
import matplotlib.pyplot as plt
import matplotlib
from datetime import datetime
from mpi4py import MPI





def read_csv_to_arrays(file_path):
    rows = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file, delimiter=";")
        rows = [row for row in csv_reader]
        return rows




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





basedir = "benchmarks/dualop_cpu_vs_gpu"
if not os.path.exists(basedir + "/graphs_amort.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)

if len(sys.argv) <= 2:
    print("not enough arguments")
    sys.exit(2)

summarize_datestr = sys.argv[1]
datestr = sys.argv[2]
graphdir = basedir + "/graphs/" + datestr
outdir = graphdir + "/amort_iters"
os.makedirs(outdir, exist_ok=True)

csv_all = read_csv_to_arrays(basedir + "/summary/" + summarize_datestr + "/timings.csv")
csv_header = csv_all[0]
csv_data = csv_all[1:]

machine_col = csv_header.index("machine")
tool_col = csv_header.index("tool")
dualop_col = csv_header.index("dual_operator")
problem_col = csv_header.index("problem")
element_type_col = csv_header.index("element_type")
n_dofs_col = csv_header.index("n_dofs")
set_col = csv_header.index("set/subdomain")
update_col = csv_header.index("update/subdomain")
apply_col = csv_header.index("apply/subdomain")
trs1storage_col = csv_header.index("TRS1_FACTOR_STORAGE")

image_index = -1













min_y_val = 1
max_y_val = 1024

problems = ["heat_transfer_2d", "heat_transfer_3d", "linear_elasticity_2d", "linear_elasticity_3d"]
# problems = ["linear_elasticity_2d", "linear_elasticity_3d"]
machines = sorted(list(set([mtd[0] for mtd in machines_tools_dualops])))
for machine in machines:
    image_index += 1
    if image_index != MPI.COMM_WORLD.Get_rank():
        continue

    tools_dualops = list(filter(lambda mtd: mtd[0] == machine, machines_tools_dualops))
    ngraphs = len(tools_dualops)

    csv_data_machine = list(filter(lambda row: row[machine_col] == machine, csv_data))

    plt.figure()
    fig, axs = plt.subplots(ngraphs, ngraphs, figsize=(ngraphs*600/100.0, ngraphs*400/100.0))

    for problem in problems:
        csv_data_problem = list(filter(lambda row: row[problem_col] == problem, csv_data_machine))
        
        physics = problem[0:-3]
        dim = int(problem[-2:-1])
        # if dim == 2: element_types = ["TRIANGLE3", "TRIANGLE6"]
        # if dim == 3: element_types = ["TETRA4", "TETRA10"]
        if dim == 2: element_types = ["TRIANGLE3"]
        if dim == 3: element_types = ["TETRA4"]
        for element_type in element_types:
            csv_data_element = list(filter(lambda row: row[element_type_col] == element_type, csv_data_problem))

            ndofs_list = sorted(list(set([int(row[n_dofs_col]) for row in csv_data_element])))

            for i in range(ngraphs):
                tool_i = tools_dualops[i][1]
                dualop_i = tools_dualops[i][2]
                csv_data_i_unsorted = list(filter(lambda row: row[tool_col] == tool_i and row[dualop_col] == dualop_i, csv_data_element))
                csv_data_i = sorted(csv_data_i_unsorted, key=lambda row: int(row[n_dofs_col]))

                for j in range(ngraphs):
                    tool_j = tools_dualops[j][1]
                    dualop_j = tools_dualops[j][2]
                    csv_data_j_unsorted = list(filter(lambda row: row[tool_col] == tool_j and row[dualop_col] == dualop_j, csv_data_element))
                    csv_data_j = sorted(csv_data_j_unsorted, key=lambda row: int(row[n_dofs_col]))

                    if i == j:
                        continue

                    x_vals = ndofs_list
                    y_vals = [float("NaN")] * len(ndofs_list)

                    for k in range(len(ndofs_list)):
                        update_i_str = csv_data_i[k][update_col]
                        update_j_str = csv_data_j[k][update_col]
                        apply_i_str = csv_data_i[k][apply_col]
                        apply_j_str = csv_data_j[k][apply_col]
                        update_i = float(update_i_str) if len(update_i_str)>0 else float("inf")
                        update_j = float(update_j_str) if len(update_j_str)>0 else float("inf")
                        apply_i = float(apply_i_str) if len(apply_i_str)>0 else float("inf")
                        apply_j = float(apply_j_str) if len(apply_j_str)>0 else float("inf")
                        is_i_update_faster = (update_i < update_j)
                        is_i_apply_faster = (apply_i < apply_j)

                        if is_i_update_faster and is_i_apply_faster:
                            # i is just faster in everything, amortization in 0 iteration
                            y = min_y_val
                        if is_i_update_faster and not is_i_apply_faster:
                            # is has faster update, but slower apply, dont plot this, its the other way around
                            y = float("NaN")
                        if not is_i_update_faster and is_i_apply_faster:
                            # i has slower update, but faster apply. There is some amortization time
                            update_diff = update_i - update_j
                            apply_diff = apply_j - apply_i
                            amortization = update_diff / apply_diff
                            y = amortization
                        if not is_i_update_faster and not is_i_apply_faster:
                            # i is just worse in everything infinite amortization time
                            y = max_y_val
                        if y == y:
                            y_vals[k] = min(max(y, min_y_val), max_y_val)
                        
                    linestyle = "-" if dim == 3 else "--"
                    color = "red" if physics == "heat_transfer" else "blue"
                    # label = str(dim) + "D-" + physics + "-" + element_type
                    label = str(dim) + "D"
                    title = tool_i + "-" + dualop_i + " VS " + tool_j + "-" + dualop_j
                    axs[i,j].loglog(x_vals, y_vals, base=2, color=color, linestyle=linestyle, label=label)
                    axs[i,j].set_title(title, fontsize="medium")
    
    xlim_min = (1 << 30)
    xlim_max = 0
    ylim_min = (1 << 30)
    ylim_max = 0
    for a in axs.flat:
        a.grid(True)
        # a.set_xlabel('n_dofs')
        # a.set_ylabel('iterations until amortization')
        if a.lines:
            a.legend(loc="upper left")
            xlim_min = min(xlim_min, a.get_xlim()[0])
            xlim_max = max(xlim_max, a.get_xlim()[1])
            ylim_min = min(ylim_min, a.get_ylim()[0])
            ylim_max = max(ylim_max, a.get_ylim()[1])
        a.set_xscale("log", base=2)
        a.set_yscale("log", base=2)
    plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
    fig.tight_layout()
    plt.savefig(outdir + "/" + machine + ".png")
    plt.close()
