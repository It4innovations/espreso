#!/usr/bin/env python3

import os
import sys
import csv
import matplotlib.pyplot as plt
import matplotlib
from datetime import datetime
from mpi4py import MPI
import mytikzplot





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
if not os.path.exists(basedir + "/graphs_time_machines.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)

if len(sys.argv) <= 2:
    print("not enough arguments")
    sys.exit(2)

summarize_datestr = sys.argv[1]
datestr = sys.argv[2]
graphdir = basedir + "/graphs/" + datestr
outdir = graphdir + "/time_machines"
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












graph_logx = 2
graph_logy = 10
graph_xlabel = "Number of DOFs per subdomain"
graph_ylabel = "Time per subdomain [ms]"
graph_legendpos = "upper left"

problems = ["heat_transfer_2d", "heat_transfer_3d", "linear_elasticity_2d", "linear_elasticity_3d"]
element_types_2d = ["TRIANGLE3", "TRIANGLE6"]
element_types_3d = ["TETRA4", "TETRA10"]
# element_types_2d = ["TRIANGLE3"]
# element_types_3d = ["TETRA4"]

machines = sorted(list(set([mtd[0] for mtd in machines_tools_dualops])))
for machine in machines:
    tools_dualops = list(filter(lambda mtd: mtd[0] == machine, machines_tools_dualops))
    
    csv_data_machine = list(filter(lambda row: row[machine_col] == machine, csv_data))
    machine_outdir = outdir + "/" + machine
    machine_tikz_outdir = outdir + "/" + machine + "/tikz"
    os.makedirs(machine_outdir, exist_ok=True)
    os.makedirs(machine_tikz_outdir, exist_ok=True)

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

            stages = ["update", "apply"]
            
            ndofs_list = sorted(list(set([int(row[n_dofs_col]) for row in csv_data_element])))
            
            plt.figure()
            fig, axs = plt.subplots(1, len(stages), figsize=(1200/100.0, 500/100.0))
            tikzplotters = mytikzplot.tikzplotter.create_table_of_plotters(1, len(stages))

            for stage_idx in range(len(stages)):
                stage = stages[stage_idx]
                if stage == "update": select_col = update_col
                if stage == "apply":  select_col = apply_col

                for mtd in tools_dualops:
                    tool = mtd[1]
                    dualop = mtd[2]

                    csv_data_tool_dualop = list(filter(lambda row: row[tool_col] == tool and row[dualop_col] == dualop, csv_data_element))
                    csv_data_sorted = sorted(csv_data_tool_dualop, key=lambda row: int(row[n_dofs_col]))

                    x_vals = ndofs_list
                    y_vals_str = [row[select_col] for row in csv_data_sorted]
                    y_vals = [float(v) if len(v)>0 else float("NaN") for v in y_vals_str]

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
                    title = machine + "-" + stage
                    tikzplotters[0][stage_idx].add_line(mytikzplot.line(x_vals, y_vals, color, linestyle, marker, label))
                    axs[stage_idx].loglog(x_vals, y_vals, base=2, color=color, linestyle=linestyle, linewidth=linewidth, marker=marker, label=label)
                    axs[stage_idx].set_title(title, fontsize="medium")

                    trs1storages = [row[trs1storage_col] for row in csv_data_sorted]
                    for i in range(len(x_vals)):
                        if trs1storages[i] == "DENSE":  axs[stage_idx].plot(x_vals[i], y_vals[i], color=color, marker="s")
                        if trs1storages[i] == "SPARSE": axs[stage_idx].plot(x_vals[i], y_vals[i], color=color, marker="P")

            xlim_min = (1 << 30)
            xlim_max = 0
            ylim_min = (1 << 30)
            ylim_max = 0
            for a in axs.flat:
                a.grid(True)
                # a.set_xlabel(graph_xlabel)
                # a.set_ylabel(graph_ylabel)
                if a.lines:
                    a.legend(loc=graph_legendpos)
                    xlim_min = min(xlim_min, a.get_xlim()[0])
                    xlim_max = max(xlim_max, a.get_xlim()[1])
                    ylim_min = min(ylim_min, a.get_ylim()[0])
                    ylim_max = max(ylim_max, a.get_ylim()[1])
                a.set_xscale("log", base=graph_logx)
                a.set_yscale("log", base=graph_logy)
            plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
            fig.tight_layout()
            plt.savefig(machine_outdir + "/" + problem + "-" + element_type + ".png")
            plt.close()

            for stage_idx in range(len(stages)):
                stage = stages[stage_idx]
                tp = tikzplotters[0][stage_idx]
                # tp.set_legendpos(graph_legendpos)
                tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
                tp.logx = graph_logx
                tp.logy = graph_logy
                tp.xlabel = graph_xlabel
                tp.ylabel = graph_ylabel
                tp.save(machine_tikz_outdir + "/" + problem + "-" + element_type + "-" + stage + ".tex")





# machines = sorted(list(set([mtd[0] for mtd in machines_tools_dualops])))
# for machine in machines:
#     tools_dualops = list(filter(lambda mtd: mtd[0] == machine, machines_tools_dualops))
    
#     csv_data_machine = list(filter(lambda row: row[machine_col] == machine, csv_data))
#     machine_outdir = outdir + "/" + machine
#     os.makedirs(machine_outdir, exist_ok=True)

#     problems = ["heat_transfer_2d", "heat_transfer_3d", "linear_elasticity_2d", "linear_elasticity_3d"]
#     for problem in problems:
#         csv_data_problem = list(filter(lambda row: row[problem_col] == problem, csv_data_machine))

#         physics = problem[0:-3]
#         dim = int(problem[-2:-1])
#         if dim == 2: element_types = ["TRIANGLE3", "TRIANGLE6"]
#         if dim == 3: element_types = ["TETRA4", "TETRA10"]
#         for element_type in element_types:
#             image_index += 1
#             if image_index % MPI.COMM_WORLD.Get_size() != MPI.COMM_WORLD.Get_rank():
#                 continue

#             csv_data_element = list(filter(lambda row: row[element_type_col] == element_type, csv_data_problem))

#             stages = ["update", "apply"]
#             approaches = ["IMPLICIT", "EXPLICIT"]
            
#             ndofs_list = sorted(list(set([int(row[n_dofs_col]) for row in csv_data_element])))
            
#             plt.figure()
#             fig, axs = plt.subplots(len(approaches), len(stages), figsize=(1920/100.0, 1080/100.0))

#             for stage_idx in range(len(stages)):
#                 stage = stages[stage_idx]
#                 if stage == "update": select_col = update_col
#                 if stage == "apply":  select_col = apply_col

#                 for approach_idx in range(len(approaches)):
#                     approach = approaches[approach_idx]

#                     csv_data_approach = list(filter(lambda row: approach in row[dualop_col], csv_data_element))
#                     tools = [mtd[1] for mtd in list(filter(lambda mtd: approach in mtd[2], tools_dualops))]

#                     for tool in tools:

#                         csv_data_tool = list(filter(lambda row: row[tool_col] == tool, csv_data_approach))
#                         csv_data_sorted = sorted(csv_data_tool, key=lambda row: int(row[n_dofs_col]))

#                         x_vals = ndofs_list
#                         y_vals_str = [row[select_col] for row in csv_data_sorted]
#                         y_vals = [float(v) if len(v)>0 else float("NaN") for v in y_vals_str]

#                         linestyle = "."
#                         if "cuda" in tool or tool == "rocm" or tool == "oneapi": linestyle = "-"
#                         else: linestyle = "--"
#                         color = "black"
#                         if tool == "cudalegacy": color = "green"
#                         if tool == "cudamodern": color = "lawngreen"
#                         if tool == "rocm": color = "darkorange"
#                         if tool == "oneapi": color = "blue"
#                         if tool == "suitesparse": color = "magenta"
#                         if tool == "mklpardiso": color = "cyan"
#                         linewidth = 2
#                         marker = None
#                         label = tool
#                         title = machine + "-" + stage + "-" + approach
#                         axs[stage_idx,approach_idx].loglog(x_vals, y_vals, base=2, color=color, linestyle=linestyle, linewidth=linewidth, marker=marker, label=label)
#                         axs[stage_idx,approach_idx].set_title(title, fontsize="medium")

#                         trs1storages = [row[trs1storage_col] for row in csv_data_sorted]
#                         for i in range(len(x_vals)):
#                             if trs1storages[i] == "DENSE":  axs[stage_idx,approach_idx].plot(x_vals[i], y_vals[i], color=color, marker="s")
#                             if trs1storages[i] == "SPARSE": axs[stage_idx,approach_idx].plot(x_vals[i], y_vals[i], color=color, marker="P")

#             xlim_min = (1 << 30)
#             xlim_max = 0
#             ylim_min = (1 << 30)
#             ylim_max = 0
#             for a in axs.flat:
#                 a.grid(True)
#                 # a.set_xlabel('n_dofs')
#                 # a.set_ylabel('time/subdomain')
#                 if a.lines:
#                     a.legend(loc="upper left")
#                     xlim_min = min(xlim_min, a.get_xlim()[0])
#                     xlim_max = max(xlim_max, a.get_xlim()[1])
#                     ylim_min = min(ylim_min, a.get_ylim()[0])
#                     ylim_max = max(ylim_max, a.get_ylim()[1])
#                 a.set_xscale("log", base=2)
#                 a.set_yscale("log", base=2)
#             plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
#             fig.tight_layout()
#             plt.savefig(machine_outdir + "/" + problem + "-" + element_type + ".png")
#             plt.close()
