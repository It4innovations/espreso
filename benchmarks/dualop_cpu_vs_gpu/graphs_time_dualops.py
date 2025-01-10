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
if not os.path.exists(basedir + "/graphs_time_dualops.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)

if len(sys.argv) <= 2:
    print("not enough arguments")
    sys.exit(2)

summarize_datestr = sys.argv[1]
datestr = sys.argv[2]
graphdir = basedir + "/graphs/" + datestr
outdir = graphdir + "/time_dualops"
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













for device in ["CPU", "GPU"]:
    problems = ["heat_transfer_2d", "heat_transfer_3d", "linear_elasticity_2d", "linear_elasticity_3d"]
    for problem in problems:
        csv_data_problem = list(filter(lambda row: row[problem_col] == problem, csv_data))

        physics = problem[0:-3]
        dim = int(problem[-2:-1])
        if dim == 2: element_types = ["TRIANGLE3", "TRIANGLE6"]
        if dim == 3: element_types = ["TETRA4", "TETRA10"]
        for element_type in element_types:
            image_index += 1
            if image_index % MPI.COMM_WORLD.Get_size() != MPI.COMM_WORLD.Get_rank():
                continue

            csv_data_element = list(filter(lambda row: row[element_type_col] == element_type, csv_data_problem))

            stages = ["update", "apply"]
            approaches = ["IMPLICIT", "EXPLICIT"]
            
            ndofs_list = sorted(list(set([int(row[n_dofs_col]) for row in csv_data_element])))
            
            plt.figure()
            fig, axs = plt.subplots(len(approaches), len(stages), figsize=(1920/100.0, 1080/100.0))

            for stage_idx in range(len(stages)):
                stage = stages[stage_idx]
                if stage == "update": select_col = update_col
                if stage == "apply":  select_col = apply_col

                for approach_idx in range(len(approaches)):
                    approach = approaches[approach_idx]

                    csv_data_approach = list(filter(lambda row: approach in row[dualop_col], csv_data_element))

                    for mtd in machines_tools_dualops:
                        machine = mtd[0]
                        tool = mtd[1]
                        dualop = mtd[2]

                        if not approach in dualop:
                            continue
                        if ("GPU" in dualop) != (device == "GPU"):
                            continue

                        csv_data_mt = list(filter(lambda row: row[machine_col] == machine and row[tool_col] == tool and row[dualop_col] == dualop, csv_data_approach))
                        csv_data_sorted = sorted(csv_data_mt, key=lambda row: int(row[n_dofs_col]))

                        x_vals = ndofs_list
                        y_vals_str = [row[select_col] for row in csv_data_sorted]
                        y_vals = [float(v) if len(v)>0 else float("NaN") for v in y_vals_str]

                        linestyle = None
                        if tool == "cudalegacy": linestyle = "-"
                        if tool == "cudamodern": linestyle = "-."
                        if tool == "rocm": linestyle = "-"
                        if tool == "oneapi": linestyle = "-"
                        if tool == "suitesparse": linestyle = "--"
                        if tool == "mklpardiso": linestyle = ":"
                        color = "black"
                        if machine == "karolina": color = "green"
                        if machine == "e4red": color = "lawngreen"
                        if machine == "lumi": color = "darkorange"
                        if machine == "tiber": color = "blue"
                        if machine == "sprddr": color = "cyan"
                        if machine == "sprhbm": color = "magenta"
                        linewidth = 2
                        marker = None
                        label = machine + "-" + tool + "-" + dualop
                        title = device + "-" + stage + "-" + approach
                        axs[stage_idx,approach_idx].loglog(x_vals, y_vals, base=2, color=color, linestyle=linestyle, linewidth=linewidth, marker=marker, label=label)
                        axs[stage_idx,approach_idx].set_title(title, fontsize="medium")

                        trs1storages = [row[trs1storage_col] for row in csv_data_sorted]
                        for i in range(len(x_vals)):
                            if trs1storages[i] == "DENSE":  axs[stage_idx,approach_idx].plot(x_vals[i], y_vals[i], color=color, marker="s")
                            if trs1storages[i] == "SPARSE": axs[stage_idx,approach_idx].plot(x_vals[i], y_vals[i], color=color, marker="P")

            xlim_min = (1 << 30)
            xlim_max = 0
            ylim_min = (1 << 30)
            ylim_max = 0
            for a in axs.flat:
                a.grid(True)
                # a.set_xlabel('n_dofs')
                # a.set_ylabel('time/subdomain [ms]')
                if a.lines:
                    a.legend(loc="upper left")
                    xlim_min = min(xlim_min, a.get_xlim()[0])
                    xlim_max = max(xlim_max, a.get_xlim()[1])
                    ylim_min = min(ylim_min, a.get_ylim()[0])
                    ylim_max = max(ylim_max, a.get_ylim()[1])
                a.set_xscale("log", base=2)
                a.set_yscale("log", base=10)
            plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
            fig.tight_layout()
            plt.savefig(outdir + "/" + device + "-" + problem + "-" + element_type + ".png")
            plt.close()
