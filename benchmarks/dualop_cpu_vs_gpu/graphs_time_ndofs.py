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
    ("karolina", "hybrid",       "EXPLICIT_SC_GPUAPPLY"),
    ("lumi",     "rocm",         "EXPLICIT_GPU"),
    ("lumi",     "rocm",         "IMPLICIT_GPU"),
    ("lumi",     "mklpardiso",   "EXPLICIT_SC"),
    ("lumi",     "mklpardiso",   "IMPLICIT"),
    ("lumi",     "suitesparse",  "EXPLICIT_SC"),
    ("lumi",     "suitesparse",  "IMPLICIT"),
    ("lumi",     "hybrid",       "EXPLICIT_SC_GPUAPPLY"),
    ("e4red",    "cudamodern",   "EXPLICIT_GPU"),
    ("e4red",    "cudamodern",   "IMPLICIT_GPU"),
    ("e4red",    "suitesparse",  "EXPLICIT_SC"),
    ("e4red",    "suitesparse",  "IMPLICIT"),
    ("e4red",    "hybrid",       "EXPLICIT_SC_GPUAPPLY"),
    ("tiber",    "oneapi",       "EXPLICIT_GPU"),
    ("tiber",    "oneapi",       "IMPLICIT_GPU"),
    ("tiber",    "mklpardiso",   "EXPLICIT_SC"),
    ("tiber",    "mklpardiso",   "IMPLICIT"),
    ("tiber",    "suitesparse",  "EXPLICIT_SC"),
    ("tiber",    "suitesparse",  "IMPLICIT"),
    ("tiber",    "hybrid",       "EXPLICIT_SC_GPUAPPLY"),
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
if not os.path.exists(basedir + "/graphs_time_ndofs.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)

if len(sys.argv) <= 2:
    print("not enough arguments")
    sys.exit(2)

summarize_datestr = sys.argv[1]
datestr = sys.argv[2]
graphdir = basedir + "/graphs/" + datestr
outdir = graphdir + "/time_ndofs"
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
min_niters = 1
max_niters = 10000

machines = sorted(list(set([mtd[0] for mtd in machines_tools_dualops])))
for machine in machines:
    tools_dualops = list(filter(lambda mtd: mtd[0] == machine, machines_tools_dualops))

    csv_data_machine = list(filter(lambda row: row[machine_col] == machine, csv_data))
    machine_outdir = outdir + "/" + machine
    os.makedirs(machine_outdir, exist_ok=True)

    for problem in problems:
        csv_data_problem = list(filter(lambda row: row[problem_col] == problem, csv_data_machine))

        physics = problem[0:-3]
        dim = int(problem[-2:-1])
        if dim == 2: element_types = element_types_2d
        if dim == 3: element_types = element_types_3d
        for element_type in element_types:
            image_index += 1
            if image_index % MPI.COMM_WORLD.Get_size() != MPI.COMM_WORLD.Get_rank():
                continue

            csv_data_element = list(filter(lambda row: row[element_type_col] == element_type, csv_data_problem))
            
            ndofs_list = sorted(list(set([int(row[n_dofs_col]) for row in csv_data_element])))
            
            n_dims = len(set([pr[-2:] for pr in problems]))
            plt.figure()
            fig, axs = plt.subplots(1, len(ndofs_list), figsize=(len(ndofs_list) * 600/100.0, 500/100.0))

            for ndofs_idx in range(len(ndofs_list)):
                ndofs = ndofs_list[ndofs_idx]

                csv_data_ndofs = list(filter(lambda row: row[n_dofs_col] == str(ndofs), csv_data_element))

                for mtd in tools_dualops:
                    tool = mtd[1]
                    dualop = mtd[2]

                    csv_data_mtd = list(filter(lambda row: row[tool_col] == tool and row[dualop_col] == dualop, csv_data_ndofs))
                    if len(csv_data_mtd) != 1:
                        print("error " + str(len(csv_data_mtd)))
                        sys.exit(1)

                    time_update_str = csv_data_mtd[0][update_col]
                    time_apply_str = csv_data_mtd[0][apply_col]
                    if len(time_update_str) == 0 or len(time_apply_str) == 0:
                        continue
                    time_update = float(time_update_str)
                    time_apply = float(time_apply_str)
                    xs = [math.exp(x / 100 * (math.log(max_niters) - math.log(min_niters)) + math.log(min_niters)) for x in range(101)]
                    ys = [time_update + x * time_apply for x in xs]

                    linestyle = ":"
                    if "IMPLICIT" in dualop: linestyle = "--"
                    if "EXPLICIT" in dualop: linestyle = "-"
                    if tool == "hybrid": linestyle = ":"
                    color = "black"
                    if tool == "cudalegacy": color = "green"
                    if tool == "cudamodern": color = "lawngreen"
                    if tool == "rocm": color = "darkorange"
                    if tool == "oneapi": color = "blue"
                    if tool == "suitesparse": color = "magenta"
                    if tool == "mklpardiso": color = "cyan"
                    if tool == "hybrid": color = "black"
                    linewidth = 2
                    marker = None
                    label = tool + "-" + dualop
                    title = "ndofs=" + str(ndofs)
                    axs[ndofs_idx].plot(xs, ys, color=color, linestyle=linestyle, linewidth=linewidth, marker=marker, label=label)
                    axs[ndofs_idx].set_title(title, fontsize="medium")

            xlim_min = (1 << 30)
            xlim_max = 0
            ylim_min = (1 << 30)
            ylim_max = 0
            for a in axs.flat:
                a.grid(True)
                # a.set_xlabel('n_dofs')
                # a.set_ylabel('time/subdomain')
                a.set_xscale("log", base=10)
                a.set_yscale("log", base=10)
                if a.lines:
                    a.legend(loc="upper left")
                    xlim_min = min(xlim_min, a.get_xlim()[0])
                    xlim_max = max(xlim_max, a.get_xlim()[1])
                    ylim_min = min(ylim_min, a.get_ylim()[0])
                    ylim_max = max(ylim_max, a.get_ylim()[1])
            plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
            fig.tight_layout()
            plt.savefig(machine_outdir + "/" + problem + "-" + element_type + ".png")
            plt.close()
