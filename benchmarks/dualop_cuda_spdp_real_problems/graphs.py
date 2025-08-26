#!/usr/bin/env python3

import os
import sys
import csv
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib
import mytikzplot



def read_csv_to_arrays(file_path):
    rows = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file, delimiter=";")
        rows = [row for row in csv_reader]
        return rows




if len(sys.argv) <= 1:
    print("not enough arguments")
    exit(1)
summary_file = sys.argv[1]

basedir = "benchmarks/dualop_cuda_spdp_real_problems"

if not os.path.isfile(summary_file):
    print("summary file/directory does not exit")
    exit(2)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")

graphs_dir = basedir + "/graphs/" + datestr
os.makedirs(graphs_dir, exist_ok=True)




csv_content = read_csv_to_arrays(summary_file)
csv_header = csv_content[0]
csv_data = csv_content[1:]

col_problem = csv_header.index("problem")
col_partstrat = csv_header.index("partition_strategy")
col_n_dofs = csv_header.index("domain_volume")
col_update_time_ps = csv_header.index("update_time_per_subdomain")
col_apply_time_ps = csv_header.index("apply_time_per_subdomain")




imgname = "speedup"
imgpath = graphs_dir + "/" + imgname + ".png"
plt.figure(figsize=(15, 10), dpi=100)

tp = mytikzplot.tikzplotter()
tppath = graphs_dir + "/" + imgname + ".tex"

problems = ["cooler", "wheel", "drillbit"]
colors = ["red", "green", "blue"]

csv_data_00 = csv_data
for i_problem in range(len(problems)):
    problem = problems[i_problem]
    csv_data_01 = [row for row in csv_data_00 if (row[col_problem] == problem)]

    csv_data_02 = sorted(csv_data_01, key=lambda row: int(row[col_n_dofs]))

    csv_data_orig = [row for row in csv_data_02 if (row[col_partstrat] == "CHUNK_COUNT")]
    csv_data_opt = [row for row in csv_data_02 if (row[col_partstrat] == "AUTO")]

    vals_x_str = [row[col_n_dofs] for row in csv_data_orig]

    vals_y_str_orig = [row[col_update_time_ps] for row in csv_data_orig]
    vals_y_str_opt = [row[col_update_time_ps] for row in csv_data_opt]

    vals_x = [float(x) for x in vals_x_str]
    vals_y_orig = [(float(y) if y != "" else float("nan")) for y in vals_y_str_orig]
    vals_y_opt = [(float(y) if y != "" else float("nan")) for y in vals_y_str_opt]

    vals_y = [float("nan")] * len(vals_x)
    for i in range(len(vals_x)):
        vals_y[i] = vals_y_orig[i] / vals_y_opt[i]

    plt.semilogx(vals_x, vals_y, base=2, color=colors[i_problem], linestyle="-", label=problem)
    tp.add_line(mytikzplot.line(vals_x, vals_y, colors[i_problem], "-", None, problem))

tp.set_bounds(plt.gca().get_xlim()[0], plt.gca().get_xlim()[1], plt.gca().get_ylim()[0], plt.gca().get_ylim()[1])
tp.logx = 2
tp.logy = 0
tp.xlabel = "Number of DOFs per subdomain"
tp.ylabel = "Speedup"
tp.save(tppath)

plt.title(imgname, fontsize="medium")
plt.legend(loc="upper left")
plt.grid(True)
plt.savefig(imgpath)
plt.close()
