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

basedir = "benchmarks/sparse_solvers_compare"

if not os.path.isfile(summary_file):
    print("summary file/directory does not exit")
    exit(2)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")

graphs_dir = basedir + "/graphs/" + datestr
os.makedirs(graphs_dir, exist_ok=True)




csv_content = read_csv_to_arrays(summary_file)
csv_header = csv_content[0]
csv_data = csv_content[1:]

col_physics = csv_header.index("physics")
col_dimension = csv_header.index("dimension")
col_element_type = csv_header.index("element_type")
col_n_dofs = csv_header.index("n_dofs")
col_update_time_ps = csv_header.index("update_time_per_subdomain")
col_apply_time_ps = csv_header.index("apply_time_per_subdomain")
col_dualop = csv_header.index("dualoperator")
col_schur = csv_header.index("schur_impl")
col_spsolver = csv_header.index("sparse_solver_impl")
col_applywhere = csv_header.index("explicit_apply_where")



tuples = [
    ("EXPLICIT_GENERALSCHUR_CPU",        "MKLPARDISO",    "AUTO",        "CPU",  "cyan", "-"),
    ("EXPLICIT_GENERALSCHUR_CPU",        "MUMPS",         "AUTO",        "CPU",  "red", "-"),
    ("EXPLICIT_GENERALSCHUR_CPU",        "SPARSE_SOLVER", "SUITESPARSE", "CPU",  "magenta", "-"),
    ("EXPLICIT_GENERALSCHUR_CPU",        "SPARSE_SOLVER", "STRUMPACK",   "CPU",  "blue", "-"),
    ("EXPLICIT_GENERALSCHUR_CPU",        "SPARSE_SOLVER", "PASTIX",      "CPU",  "green", "-"),
    ("EXPLICIT_GENERALSCHUR_CPU",        "TRIANGULAR",    "AUTO",        "CPU",  "black", "-."),
    ("EXPLICIT_GENERALSCHUR_CPU",        "AUTO",          "AUTO",        "GPU",  "black", ":"),
    ("IMPLICIT_GENERALSPARSESOLVER_CPU", "AUTO",          "MKLPARDISO",  "AUTO", "cyan", "--"),
    ("IMPLICIT_GENERALSPARSESOLVER_CPU", "AUTO",          "SUITESPARSE", "AUTO", "magenta", "--"),
    ("IMPLICIT_GENERALSPARSESOLVER_CPU", "AUTO",          "MUMPS",       "AUTO", "red", "--"),
    ("IMPLICIT_GENERALSPARSESOLVER_CPU", "AUTO",          "STRUMPACK",   "AUTO", "blue", "--"),
    ("IMPLICIT_GENERALSPARSESOLVER_CPU", "AUTO",          "PASTIX",      "AUTO", "green", "--")
]

physics_list = ["heat_transfer"]
dimensions_list = ["2", "3"]
element_types_2d_list = ["TRIANGLE3"]
element_types_3d_list = ["TETRA4"]

subplots_counts = [2,1]
my_figsize_x = 4000
my_figsize_y = 3000

image_index = -1

csv_data_00 = csv_data
for physics in physics_list:
    csv_data_02 = [row for row in csv_data_00 if (row[col_physics] == physics)]
    for dimension in dimensions_list:
        csv_data_03 = [row for row in csv_data_02 if (row[col_dimension] == dimension)]
        if dimension == "2": element_types_list = element_types_2d_list
        if dimension == "3": element_types_list = element_types_3d_list
        for element_type in element_types_list:
            csv_data_04 = [row for row in csv_data_03 if (row[col_element_type] == element_type)]

            upd_app_list = ["update", "apply"]
            for ua in upd_app_list:
                if ua == "update":
                    col_data = col_update_time_ps
                if ua == "apply":
                    col_data = col_apply_time_ps

                imgname = physics + "-" + dimension + "D-" + element_type + "-" + ua
                imgpath = graphs_dir + "/" + imgname + ".png"
                plt.figure(figsize=(15, 10), dpi=100)

                tp = mytikzplot.tikzplotter()
                tppath = graphs_dir + "/" + imgname + ".tikz"

                for tupleval in tuples:
                    dualop = tupleval[0]
                    schur = tupleval[1]
                    spsolver = tupleval[2]
                    applywhere = tupleval[3]
                    line_color = tupleval[4]
                    line_style = tupleval[5]

                    tuple_name = dualop + "-" + schur + "-" + spsolver + "-" + applywhere

                    csv_data_05 = [row for row in csv_data_04 if (row[col_dualop] == dualop)]
                    csv_data_06 = [row for row in csv_data_05 if (row[col_schur] == schur)]
                    csv_data_07 = [row for row in csv_data_06 if (row[col_spsolver] == spsolver)]
                    csv_data_08 = [row for row in csv_data_07 if (row[col_applywhere] == applywhere)]
                    csv_data_09 = sorted(csv_data_08, key=lambda row: int(row[col_n_dofs]))
        
                    vals_x_str = [row[col_n_dofs] for row in csv_data_09]
                    vals_y_str = [row[col_data] for row in csv_data_09]
                    vals_x = [float(x) for x in vals_x_str]
                    vals_y = [(float(y) if y != "" else float("nan")) for y in vals_y_str]

                    plt.loglog(vals_x, vals_y, base=2, color=line_color, linestyle=line_style, label=tuple_name)
                    tp.add_line(mytikzplot.line(vals_x, vals_y, line_color, line_style, None, tuple_name))

                plt.title(imgname, fontsize="medium")
                plt.legend(loc="upper left")
                plt.grid(True)
                plt.savefig(imgpath)
                plt.close()

                tp.save(tppath)

