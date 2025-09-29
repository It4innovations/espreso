#!/usr/bin/env python3

import os
import sys
import csv
import math
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

basedir = "benchmarks/cluster_domain_size_2"

if not os.path.isfile(summary_file):
    print("summary file/directory does not exit")
    exit(2)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")

graphs_dir = basedir + "/graphs/" + datestr
graphs_dir_png = graphs_dir + "/png_speedup_2"
graphs_dir_tikz = graphs_dir + "/tikz_speedup_2"
os.makedirs(graphs_dir_png, exist_ok=True)
os.makedirs(graphs_dir_tikz, exist_ok=True)




csv_content = read_csv_to_arrays(summary_file)
csv_header = csv_content[0]
csv_data = csv_content[1:]

col_physics = csv_header.index("physics")
col_dimension = csv_header.index("dimension")
col_element_type = csv_header.index("element_type")
col_dualop = csv_header.index("dual_operator")
col_precond = csv_header.index("preconditioner")
col_approach = csv_header.index("approach")
col_n_dofs_cluster = csv_header.index("n_dofs_total_calc")
col_n_domains_per_cluster = csv_header.index("domains_per_cluster")
col_n_dofs_per_domain = csv_header.index("n_dofs_domain_calc")
col_dualop_update_time = csv_header.index("dualop_update_time")
col_dualop_apply_time = csv_header.index("dualop_apply_time")
col_precond_update_time = csv_header.index("precond_update_time")
col_precond_apply_time = csv_header.index("precond_apply_time")
col_total_update_time = csv_header.index("total_update_time")
col_total_apply_time = csv_header.index("total_apply_time")
col_feti_method = csv_header.index("feti_method")



# physics_list = sorted(list(set([row[col_physics] for row in csv_data])))
# element_types_list = sorted(list(set([row[col_element_type] for row in csv_data])))
# approach_list = sorted(list(set([row[col_approach] for row in csv_data])))
# feti_method_list = sorted(list(set([row[col_feti_method] for row in csv_data])))
physics_list = ["heat_transfer", "linear_elasticity"]
element_types_list = ["TRIANGLE3", "TETRA4"]
approach_list = ["IMPLICIT_CPU", "EXPLICIT_CPU", "EXPLICIT_GPU"]
feti_method_list = ["TOTAL_FETI", "HYBRID_FETI"]

division_variant_cols = [col_n_domains_per_cluster, col_n_dofs_per_domain]
division_variant_labels = ["Number of subdomains per cluster", "Number of DOFs per subdomain"]
division_variant_labels_lite = ["domain_count", "domain_size"]

part_labels = ["dualop", "precond", "total"]
phase_labels = ["update", "apply"]
part_phase_cols = [[col_dualop_update_time, col_dualop_apply_time], [col_precond_update_time, col_precond_apply_time], [col_total_update_time, col_total_apply_time]]

my_figsize_x = 500
my_figsize_y = 400

image_index = -1

niters_min = 1
niters_max = 10000
niters_intervals = 100
niters_list = range(niters_intervals+1)
niters_list = [x / niters_intervals for x in niters_list]
niters_list = [x * (math.log(niters_max) - math.log(niters_min)) for x in niters_list]
niters_list = [x + math.log(niters_min) for x in niters_list]
niters_list = [math.exp(x) for x in niters_list]

csv_data_00 = csv_data
for physics in physics_list:
    for element_type in element_types_list:
        for division_variant_idx in range(len(division_variant_cols)):
            division_variant_col = division_variant_cols[division_variant_idx]
            division_variant = division_variant_labels[division_variant_idx]
            division_variant_lite = division_variant_labels_lite[division_variant_idx]
            image_index = image_index + 1
            if (image_index % MPI.COMM_WORLD.Get_size()) != MPI.COMM_WORLD.Get_rank():
                continue

            csv_data_02 = [row for row in csv_data_00 if ((row[col_physics] == physics) and (row[col_element_type] == element_type))]
            csv_data_03 = sorted(csv_data_02, key=lambda row: int(row[division_variant_col]))
            dimension = csv_data_03[0][col_dimension]

            n_dofs_cluster_list = sorted(list(set([row[col_n_dofs_cluster] for row in csv_data_03])), key=lambda x: int(x))
            grafs_nrows = len(part_labels) * len(feti_method_list)
            grafs_ncols = len(n_dofs_cluster_list)

            imgname = physics + "-" + dimension + "D-" + element_type + "-" + division_variant_lite
            imgpath = graphs_dir_png + "/" + imgname + ".png"
            plt.figure()
            fig, axs = plt.subplots(grafs_nrows, grafs_ncols, figsize=(my_figsize_x*grafs_ncols/100.0, my_figsize_y*grafs_nrows/100.0), sharex=True, sharey=True)

            tps = mytikzplot.tikzplotter.create_table_of_plotters(grafs_nrows, grafs_ncols)

            for n_dofs_cluster_idx in range(len(n_dofs_cluster_list)):
                n_dofs_cluster = n_dofs_cluster_list[n_dofs_cluster_idx]
                csv_data_04 = [row for row in csv_data_03 if (row[col_n_dofs_cluster] == n_dofs_cluster)]

                for feti_method_idx in range(len(feti_method_list)):
                    feti_method = feti_method_list[feti_method_idx]
                    csv_data_06 = [row for row in csv_data_04 if (row[col_feti_method] == feti_method)]

                    for part_idx in range(len(part_labels)):
                        part_label = part_labels[part_idx]

                        graph_row = part_idx + len(part_labels) * (feti_method_idx)
                        graph_col = n_dofs_cluster_idx

                        title = part_label + ", " + feti_method + ", " + n_dofs_cluster + " DOFs in cluster"
                        axs[graph_row,graph_col].set_title(title)
                        tps[graph_row][graph_col].title = title

                        division_vals_raw = list(set([row[division_variant_col] for row in csv_data_06]))
                        division_vals = sorted(division_vals_raw, key=lambda x: int(x))

                        for division_val_idx in range(len(division_vals)):
                            division_val = division_vals[division_val_idx]
                            csv_data_07 = [row for row in csv_data_06 if (row[division_variant_col] == division_val)]

                            time_ic_update_str = [row[part_phase_cols[part_idx][0]] for row in csv_data_07 if (row[col_approach] == "IMPLICIT_CPU")]
                            time_eg_update_str = [row[part_phase_cols[part_idx][0]] for row in csv_data_07 if (row[col_approach] == "EXPLICIT_GPU")]
                            time_ic_apply_str = [row[part_phase_cols[part_idx][1]] for row in csv_data_07 if (row[col_approach] == "IMPLICIT_CPU")]
                            time_eg_apply_str = [row[part_phase_cols[part_idx][1]] for row in csv_data_07 if (row[col_approach] == "EXPLICIT_GPU")]
                            if len(time_ic_update_str) != 1 or len(time_eg_update_str) != 1 or len(time_ic_apply_str) != 1 or len(time_eg_apply_str) != 1:
                                print("something is wrong")
                            time_ic_update_str = time_ic_update_str[0]
                            time_eg_update_str = time_eg_update_str[0]
                            time_ic_apply_str = time_ic_apply_str[0]
                            time_eg_apply_str = time_eg_apply_str[0]
                            time_ic_update = float("nan") if len(time_ic_update_str) == 0 else float(time_ic_update_str)
                            time_eg_update = float("nan") if len(time_eg_update_str) == 0 else float(time_eg_update_str)
                            time_ic_apply = float("nan") if len(time_ic_apply_str) == 0 else float(time_ic_apply_str)
                            time_eg_apply = float("nan") if len(time_eg_apply_str) == 0 else float(time_eg_apply_str)

                            label = division_variant_lite + "=" + str(division_val)
                            vals_x = niters_list
                            vals_y = [float("nan")] * len(niters_list)

                            for i in range(len(niters_list)):
                                niters = niters_list[i]
                                time_ic = time_ic_update + niters * time_ic_apply
                                time_eg = time_eg_update + niters * time_eg_apply
                                vals_y[i] = time_ic / time_eg

                            color = ["red", "cyan", "green", "magenta", "blue"][division_val_idx % 5]
                            linestyle = ["-", "--", ":"][division_val_idx % 3]
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
                    a.legend(loc="upper left", fontsize="x-small")
                    a.set_xlabel("Number of applications")
                    a.set_ylabel("Speedup")
                    tp.xlabel = "Number of applications"
                    tp.ylabel = "Speedup"
                    [xlim_min, xlim_max] = a.get_xlim()
                    [ylim_min, ylim_max] = a.get_ylim()
                    tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
            fig.tight_layout()
            plt.savefig(imgpath)
            plt.close()

            for graph_row in range(grafs_nrows):
                for graph_col in range(grafs_ncols):
                    path = graphs_dir_tikz + "/" + imgname + "-" + part_labels[graph_row % 3] + "-" + feti_method_list[graph_row // 3] + "-" + n_dofs_cluster_list[graph_col] + ".tex"
                    tps[graph_row][graph_col].save(path)
