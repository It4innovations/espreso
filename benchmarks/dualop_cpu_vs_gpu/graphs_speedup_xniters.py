#!/usr/bin/env python3

import os
import sys
import csv
import math
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
if not os.path.exists(basedir + "/graphs_speedup_xniters.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)

if len(sys.argv) <= 2:
    print("not enough arguments")
    sys.exit(2)

summarize_datestr = sys.argv[1]
datestr = sys.argv[2]
graphdir = basedir + "/graphs/" + datestr
outdir = graphdir + "/speedup_xniters"
tikzoutdir = outdir + "/tikz"
os.makedirs(outdir, exist_ok=True)
os.makedirs(tikzoutdir, exist_ok=True)

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

















graph_logx = 10
graph_logy = 10
graph_xlabel = "Number of iterations"
graph_ylabel = "Speedup"
graph_legendpos = "upper left"

problems = ["heat_transfer_2d", "heat_transfer_3d", "linear_elasticity_2d", "linear_elasticity_3d"]
# element_types_2d = ["TRIANGLE3", "TRIANGLE6"]
# element_types_3d = ["TETRA4", "TETRA10"]
element_types_2d = ["TRIANGLE3"]
element_types_3d = ["TETRA4"]
min_niters = 1
max_niters = 10000
step = 1.1
nintervals = int((math.log(max_niters) - math.log(min_niters)) / math.log(step)) + 1
npoints = nintervals + 1
niters_list = [math.exp(x / nintervals * (math.log(max_niters) - math.log(min_niters)) + math.log(min_niters)) for x in range(npoints)]

machines = sorted(list(set([mtd[0] for mtd in machines_tools_dualops])))
for machine in machines:
    image_index += 1
    if image_index % MPI.COMM_WORLD.Get_size() != MPI.COMM_WORLD.Get_rank():
        continue

    csv_data_machine = list(filter(lambda row: row[machine_col] == machine, csv_data))

    tools_dualops = list(filter(lambda mtd: mtd[0] == machine, machines_tools_dualops))
            
    n_dims = len(set([pr[-2:] for pr in problems]))
    plt.figure()
    ngraphs_x = len(element_types_2d) * len(problems) // n_dims
    ngraphs_y = n_dims
    fig, axs = plt.subplots(ngraphs_y, ngraphs_x, figsize=(ngraphs_x * 700/100.0, n_dims * 500/100.0))
    tikzplotters = mytikzplot.tikzplotter.create_table_of_plotters(ngraphs_y, ngraphs_x)

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
            color = "red"

            for ndofs_idx in range(len(ndofs_list)):
                ndofs = ndofs_list[ndofs_idx]
                csv_data_ndofs = list(filter(lambda row: row[n_dofs_col] == str(ndofs), csv_data_element))
                ys = [None] * len(niters_list)

                for niters_idx in range(len(niters_list)):
                    niters = niters_list[niters_idx]
                    best_mtd = None
                    best_time = 1000000000

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
                        time = time_update + niters * time_apply
                        if time < best_time:
                            best_time = time
                            best_mtd = mtd
                    
                    orig_mtd = (machine, "mklpardiso", "IMPLICIT")
                    orig_time = -1
                    if not orig_mtd in machines_tools_dualops:
                        orig_mtd = (machine, "suitesparse", "IMPLICIT")
                    if True:
                        csv_data_mtd = list(filter(lambda row: row[tool_col] == orig_mtd[1] and row[dualop_col] == orig_mtd[2], csv_data_ndofs))
                        if len(csv_data_mtd) != 1:
                            print("error " + str(len(csv_data_mtd)))
                            sys.exit(1)
                        time_update_str = csv_data_mtd[0][update_col]
                        time_apply_str = csv_data_mtd[0][apply_col]
                        if len(time_update_str) > 0 and len(time_apply_str) > 0:
                            time_update = float(time_update_str)
                            time_apply = float(time_apply_str)
                            orig_time = time_update + niters * time_apply

                    if orig_time == -1:
                        speedup = 1000
                    else:
                        speedup = orig_time / best_time
                    ys[niters_idx] = speedup

                if "2d" in problem: axs_y = 0
                else: axs_y = 1
                axs_x = element_type_idx
                if "linear_elasticity" in problem: axs_x += len(element_types)
                if color == "red": color = "blue"
                elif color == "blue": color = "green"
                else: color = "red"
                linestyle = "-"
                linewidth = 2
                title = problem + "-" + element_type
                xs = list(niters_list)
                last_one = (len(ys) - ys[::-1].index(1.0) - 1) if (1.0 in ys) else 0
                xs = xs[last_one:]
                ys = ys[last_one:]
                tikzplotters[axs_y][axs_x].add_line(mytikzplot.line(xs, ys, color, linestyle, None, None))
                axs[axs_y,axs_x].plot(xs, ys, color=color, linestyle=linestyle, linewidth=linewidth)
                axs[axs_y,axs_x].set_title(title, fontsize="medium")
                if niters == niters_list[-1]:
                    text = str(ndofs) + " DOFs"
                    x = niters_list[-1] * 1.1
                    y = ys[-1]
                    axs[axs_y,axs_x].text(x, y, text, fontsize="medium", verticalalignment="center", bbox=dict(edgecolor="none", facecolor="white", alpha=0.7, pad=0))
                    tikzplotters[axs_y][axs_x].add_text(mytikzplot.textbox(text, x, y))

    xlim_min = (1 << 30)
    xlim_max = 0
    ylim_min = (1 << 30)
    ylim_max = 0
    for a in axs.flat:
        a.grid(True)
        # a.set_xlabel(graph_xlabel)
        # a.set_ylabel(graph_ylabel)
        a.set_xscale("log", base=graph_logx)
        a.set_yscale("log", base=graph_logy)
        if a.lines:
            # a.legend(loc=graph_legendpos)
            xlim_min = min(xlim_min, a.get_xlim()[0])
            xlim_max = max(xlim_max, a.get_xlim()[1])
            ylim_min = min(ylim_min, a.get_ylim()[0])
            ylim_max = max(ylim_max, a.get_ylim()[1])
    xlim_max *= 5
    plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
    fig.tight_layout()
    plt.savefig(outdir + "/" + machine + ".png")
    plt.close()

    for problem in problems:
        physics = problem[0:-3]
        dim = int(problem[-2:-1])
        if dim == 2: element_types = element_types_2d
        if dim == 3: element_types = element_types_3d
        for element_type_idx in range(len(element_types)):
            element_type = element_types[element_type_idx]
            if "2d" in problem: axs_y = 0
            else: axs_y = 1
            axs_x = element_type_idx
            if "linear_elasticity" in problem: axs_x += len(element_types)
            tp = tikzplotters[axs_y][axs_x]
            # tp.set_legendpos(graph_legendpos)
            tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
            tp.logx = graph_logx
            tp.logy = graph_logy
            tp.xlabel = graph_xlabel
            tp.ylabel = graph_ylabel
            tp.save(tikzoutdir + "/" + machine + "-" + problem + "-" + element_type + ".tex")
