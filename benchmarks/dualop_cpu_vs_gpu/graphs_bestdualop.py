#!/usr/bin/env python3

import os
import sys
import csv
import math
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
if not os.path.exists(basedir + "/graphs_bestdualop.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)

if len(sys.argv) <= 2:
    print("not enough arguments")
    sys.exit(2)

summarize_datestr = sys.argv[1]
datestr = sys.argv[2]
graphdir = basedir + "/graphs/" + datestr
outdir = graphdir + "/bestdualop"
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













problems = ["heat_transfer_2d", "heat_transfer_3d", "linear_elasticity_2d", "linear_elasticity_3d"]
element_types_2d = ["TRIANGLE3", "TRIANGLE6"]
element_types_3d = ["TETRA4", "TETRA10"]
# element_types_2d = ["TRIANGLE3"]
# element_types_3d = ["TETRA4"]
max_niters = 10000

machines = sorted(list(set([mtd[0] for mtd in machines_tools_dualops])))
for machine in machines:
    image_index += 1
    if image_index != MPI.COMM_WORLD.Get_rank():
        continue

    csv_data_machine = list(filter(lambda row: row[machine_col] == machine, csv_data))

    tools_dualops = list(filter(lambda mtd: mtd[0] == machine, machines_tools_dualops))
            
    n_dims = len(set([pr[-2:] for pr in problems]))
    plt.figure()
    ngraphs_x = len(element_types_2d) * len(problems) // n_dims
    fig, axs = plt.subplots(n_dims, ngraphs_x, figsize=(ngraphs_x * 700/100.0, n_dims * 500/100.0))

    for problem in problems:
        csv_data_problem = list(filter(lambda row: row[problem_col] == problem, csv_data_machine))

        physics = problem[0:-3]
        dim = int(problem[-2:-1])
        if dim == 2: element_types = element_types_2d
        if dim == 3: element_types = element_types_3d
        for element_type_idx in range(len(element_types)):
            element_type = element_types[element_type_idx]
            csv_data_element = list(filter(lambda row: row[element_type_col] == element_type, csv_data_problem))
            
            ndofs_list = sorted(list(set([int(row[n_dofs_col]) for row in csv_data_element])))
            prev_ndofs_label = 0

            already_plotted_mtds = []

            for ndofs in ndofs_list:
                csv_data_ndofs = list(filter(lambda row: row[n_dofs_col] == str(ndofs), csv_data_element))

                curr_niters = 1
                curr_time = 1000000000
                curr_update = 1000000000
                curr_apply = 1000000000
                curr_mtd = 0

                for mtd in tools_dualops:
                    csv_data_mtd = list(filter(lambda row: row[tool_col] == mtd[1] and row[dualop_col] == mtd[2], csv_data_ndofs))
                    if len(csv_data_mtd) != 1:
                        print("error " + str(len(csv_data_mtd)))
                        sys.exit(1)
                    time_update_str = csv_data_mtd[0][update_col]
                    time_apply_str = csv_data_mtd[0][apply_col]
                    if len(time_update_str) == 0 or len(time_apply_str) == 0: continue
                    time_update = float(time_update_str)
                    time_apply = float(time_apply_str)
                    time = time_update + curr_niters * time_apply
                    if time < curr_time:
                        curr_time = time
                        curr_mtd = mtd
                        curr_update = time_update
                        curr_apply = time_apply
                
                while True:
                    next_niters = 1000000000
                    next_mtd = 0
                    next_update = 1000000000
                    next_apply = 1000000000
                    next_time = 1000000000
                    for mtd in tools_dualops:
                        if mtd == curr_mtd:
                            continue
                        csv_data_mtd = list(filter(lambda row: row[tool_col] == mtd[1] and row[dualop_col] == mtd[2], csv_data_ndofs))
                        time_update_str = csv_data_mtd[0][update_col]
                        time_apply_str = csv_data_mtd[0][apply_col]
                        if len(time_update_str) == 0 or len(time_apply_str) == 0:
                            continue
                        time_update = float(time_update_str)
                        time_apply = float(time_apply_str)
                        if time_apply >= curr_apply:
                            continue
                        niters = (time_update - curr_update) / (curr_apply - time_apply)
                        if niters > curr_niters and niters < next_niters:
                            next_niters = niters
                            next_mtd = mtd
                            next_update = time_update
                            next_apply = time_apply
                    if next_niters >= max_niters:
                        next_niters = max_niters
                    next_time = curr_update + next_niters * curr_apply

                    if "2d" in problem: axs_y = 0
                    else: axs_y = 1
                    axs_x = element_type_idx
                    if "linear_elasticity" in problem: axs_x += len(element_types)
                    color = "black"
                    if curr_mtd[1] == "cudalegacy": color = "green"
                    if curr_mtd[1] == "cudamodern": color = "lawngreen"
                    if curr_mtd[1] == "rocm": color = "darkorange"
                    if curr_mtd[1] == "oneapi": color = "blue"
                    if curr_mtd[1] == "suitesparse": color = "magenta"
                    if curr_mtd[1] == "mklpardiso": color = "cyan"
                    linestyle = ":"
                    if "IMPLICIT" in curr_mtd[2]: linestyle = "--"
                    if "EXPLICIT" in curr_mtd[2]: linestyle = "-"
                    linewidth = 2
                    marker = None
                    label = curr_mtd[1] + "-" + curr_mtd[2]
                    if label in already_plotted_mtds: label = None
                    else: already_plotted_mtds.append(label)
                    title = problem + "-" + element_type
                    xs = [math.exp(x / 100 * (math.log(next_niters) - math.log(curr_niters)) + math.log(curr_niters)) for x in range(101)]
                    ys = [x * curr_apply + curr_update for x in xs]
                    axs[axs_y,axs_x].plot(xs, ys, color=color, linestyle=linestyle, linewidth=linewidth, marker=marker, label=label)
                    axs[axs_y,axs_x].set_title(title, fontsize="medium")
                    if next_niters >= max_niters and ys[-1] > 1.7 * prev_ndofs_label:
                        axs[axs_y,axs_x].text(xs[-1] * 1.1, ys[-1], str(ndofs) + " DOFs", fontsize="medium", verticalalignment="center", bbox=dict(edgecolor="none", facecolor="white", alpha=0.7, pad=0))
                        prev_ndofs_label = ys[-1]

                    if next_niters >= max_niters:
                        break
                    curr_niters = next_niters
                    curr_update = next_update
                    curr_apply = next_apply
                    curr_time = next_time
                    curr_mtd = next_mtd

    xlim_min = (1 << 30)
    xlim_max = 0
    ylim_min = (1 << 30)
    ylim_max = 0
    for a in axs.flat:
        a.grid(True)
        # a.set_xlabel('n_iters')
        # a.set_ylabel('time/subdomain')
        a.set_xscale("log", base=10)
        a.set_yscale("log", base=10)
        if a.lines:
            a.legend(loc="upper left")
            xlim_min = min(xlim_min, a.get_xlim()[0])
            xlim_max = max(xlim_max, a.get_xlim()[1])
            ylim_min = min(ylim_min, a.get_ylim()[0])
            ylim_max = max(ylim_max, a.get_ylim()[1])
    xlim_max *= 5
    plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
    fig.tight_layout()
    plt.savefig(outdir + "/" + machine + ".png")
    plt.close()
