#!/usr/bin/env python3

import os
import sys
import csv
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib
import mytikzplot
from mpi4py import MPI



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

basedir = "benchmarks/dualop_cluster_domain_size"

if not os.path.isfile(summary_file):
    print("summary file/directory does not exit")
    exit(2)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")

graphs_dir = basedir + "/graphs/" + datestr
graphs_dir_png = graphs_dir + "/png_times_combined2"
graphs_dir_tikz = graphs_dir + "/tikz_times_combined2"
os.makedirs(graphs_dir_png, exist_ok=True)
os.makedirs(graphs_dir_tikz, exist_ok=True)




csv_content = read_csv_to_arrays(summary_file)
csv_header = csv_content[0]
csv_data = csv_content[1:]

col_physics = csv_header.index("physics")
col_dimension = csv_header.index("dimension")
col_element_type = csv_header.index("element_type")
col_dualop = csv_header.index("dual_operator")
col_n_dofs_cluster = csv_header.index("n_dofs_total_calc")
col_n_domains_per_cluster = csv_header.index("domains_per_cluster")
col_n_dofs_per_domain = csv_header.index("n_dofs_domain_calc")
col_update_time = csv_header.index("update_time")
col_apply_time = csv_header.index("apply_time")
col_feti_method = csv_header.index("feti_method")



# physics_list = sorted(list(set([row[col_physics] for row in csv_data])))
# element_types_list = sorted(list(set([row[col_element_type] for row in csv_data])))
# dual_operator_list = sorted(list(set([row[col_dualop] for row in csv_data])))
# feti_method_list = sorted(list(set([row[col_feti_method] for row in csv_data])))
physics_list = ["heat_transfer", "linear_elasticity"]
element_types_list = ["TRIANGLE3", "TRIANGLE6", "SQUARE4", "SQUARE8", "TETRA4", "TETRA10", "HEXA8", "HEXA20"]
dual_operator_list = ["IMPLICIT_GENERALSPARSESOLVER_CPU", "EXPLICIT_GENERALSCHUR_CPU", "EXPLICIT_GENERALSCHUR_GPU"]
feti_method_list = ["TOTAL_FETI", "HYBRID_FETI"]


subdivision_variant_cols = [col_n_domains_per_cluster, col_n_dofs_per_domain]
subdivision_variant_labels = ["Number of subdomains per cluster", "Number of DOFs per subdomain"]

update_apply_variant_cols = [col_update_time, col_apply_time]
update_apply_variant_labels = ["preprocessing", "application"]

my_figsize_x = 500
my_figsize_y = 400

image_index = -1

niters_list = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000]
niters_list_linecolor = ["blue", "red", "green", "orange", "blue", "red", "green", "orange", "blue"]
niters_list_linestyle = ["-", "--", ":", "-", "--", ":", "-", "--", ":"]

csv_data_00 = csv_data
for physics in physics_list:
    for element_type in element_types_list:
        for dualop in dual_operator_list:
            image_index = image_index + 1
            if (image_index % MPI.COMM_WORLD.Get_size()) != MPI.COMM_WORLD.Get_rank():
                continue

            csv_data_03 = [row for row in csv_data_00 if ((row[col_physics] == physics) and (row[col_element_type] == element_type) and (row[col_dualop] == dualop))]
            dimension = csv_data_03[0][col_dimension]

            n_dofs_cluster_list = sorted(list(set([row[col_n_dofs_cluster] for row in csv_data_03])), key=lambda x: int(x))
            grafs_nrows = len(subdivision_variant_cols) * len(feti_method_list)
            grafs_ncols = len(n_dofs_cluster_list)

            imgname = physics + "-" + dimension + "D-" + element_type + "-" + dualop
            imgpath = graphs_dir_png + "/" + imgname + ".png"
            plt.figure()
            fig, axs = plt.subplots(grafs_nrows, grafs_ncols, figsize=(my_figsize_x*grafs_ncols/100.0, my_figsize_y*grafs_nrows/100.0), sharex="row", sharey=True)

            tps = mytikzplot.tikzplotter.create_table_of_plotters(grafs_nrows, grafs_ncols)

            for n_dofs_cluster_idx in range(len(n_dofs_cluster_list)):
                n_dofs_cluster = n_dofs_cluster_list[n_dofs_cluster_idx]
                csv_data_04 = [row for row in csv_data_03 if (row[col_n_dofs_cluster] == n_dofs_cluster)]

                for subdivision_variant_idx in range(len(subdivision_variant_cols)):
                    subdivision_variant_col = subdivision_variant_cols[subdivision_variant_idx]
                    csv_data_05 = sorted(csv_data_04, key=lambda row: int(row[subdivision_variant_col]))

                    for feti_method_idx in range(len(feti_method_list)):
                        feti_method = feti_method_list[feti_method_idx]
                        csv_data_06 = [row for row in csv_data_05 if (row[col_feti_method] == feti_method)]

                        graph_row = 2 * subdivision_variant_idx + feti_method_idx
                        graph_col = n_dofs_cluster_idx

                        axs[graph_row,graph_col].set_title(feti_method + ", " + n_dofs_cluster + " DOFs in cluster")

                        vals_x_str = [row[subdivision_variant_col] for row in csv_data_06]
                        vals_x = [float(x) for x in vals_x_str]

                        update_times_str = [row[col_update_time] for row in csv_data_06]
                        apply_times_str = [row[col_apply_time] for row in csv_data_06]
                        update_times = [float("nan") if len(v)==0 else float(v) for v in update_times_str]
                        apply_times = [float("nan") if len(v)==0 else float(v) for v in apply_times_str]

                        for niters_idx in range(len(niters_list)):
                            niters = niters_list[niters_idx]

                            vals_y = [float("nan")] * len(vals_x)
                            for i in range(len(vals_y)):
                                vals_y[i] = update_times[i] + niters * apply_times[i]

                            color = niters_list_linecolor[niters_idx]
                            linestyle = niters_list_linestyle[niters_idx]
                            label = str(niters) + " iterations"

                            axs[graph_row,graph_col].plot(vals_x, vals_y, color=color, linestyle=linestyle, label=label)
                            tps[graph_row][graph_col].add_line(mytikzplot.line(vals_x, vals_y, color, linestyle, None, label))



            for graph_row in range(grafs_nrows):
                for graph_col in range(grafs_ncols):
                    tp = tps[graph_row][graph_col]
                    a = axs[graph_row, graph_col]
                    a.grid(True)
                    a.xaxis.set_tick_params(labelbottom=True)
                    a.yaxis.set_tick_params(labelleft=True)
                    a.set_xscale("log", base=10)
                    a.set_yscale("log", base=10)
                    tp.logx = 10
                    tp.logy = 10
                    if graph_row == 0 and graph_col == 0:
                        a.legend(loc="upper left")
                    a.set_xlabel(subdivision_variant_labels[graph_row // 2])
                    a.set_ylabel("Step time [ms]")
                    tp.xlabel = subdivision_variant_labels[graph_row // 2]
                    tp.ylabel = "Step time [ms]"
                    [xlim_min, xlim_max] = a.get_xlim()
                    [ylim_min, ylim_max] = a.get_ylim()
                    tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
            fig.tight_layout()
            plt.savefig(imgpath)
            plt.close()

            for graph_row in range(grafs_nrows):
                for graph_col in range(grafs_ncols):
                    path = graphs_dir_tikz + "/" + imgname + "-" + feti_method_list[graph_row % 2] + "-" + subdivision_variant_labels[graph_row // 2] + "-" + n_dofs_cluster_list[graph_col] + ".tex"
                    tps[graph_row][graph_col].save(path)
