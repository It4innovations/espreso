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

basedir = "benchmarks/cluster_domain_size_2"

if not os.path.isfile(summary_file):
    print("summary file/directory does not exit")
    exit(2)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")

graphs_dir = basedir + "/graphs/" + datestr
graphs_dir_png = graphs_dir + "/png_speedup_1"
graphs_dir_tikz = graphs_dir + "/tikz_speedup_1"
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

                        vals_x_str = [row[division_variant_col] for row in csv_data_06 if row[col_approach] == "IMPLICIT_CPU"]
                        vals_x = [float(x) for x in vals_x_str]

                        for phase_idx in range(len(phase_labels)):
                            phase_label = phase_labels[phase_idx]
                            value_col = part_phase_cols[part_idx][phase_idx]
                            
                            times_ic_str = [row[value_col] for row in csv_data_06 if (row[col_approach] == "IMPLICIT_CPU")]
                            times_ec_str = [row[value_col] for row in csv_data_06 if (row[col_approach] == "EXPLICIT_CPU")]
                            times_eg_str = [row[value_col] for row in csv_data_06 if (row[col_approach] == "EXPLICIT_GPU")]
                            times_ic = [float(val) if len(val) > 0 else float("nan") for val in times_ic_str]
                            times_ec = [float(val) if len(val) > 0 else float("nan") for val in times_ec_str]
                            times_eg = [float(val) if len(val) > 0 else float("nan") for val in times_eg_str]

                            spdn_eg_ic = [float("nan")] * len(vals_x)
                            spdn_eg_ec = [float("nan")] * len(vals_x)
                            spdn_ec_ic = [float("nan")] * len(vals_x)
                            for i in range(len(vals_x)):
                                spdn_eg_ic[i] = times_ic[i] / times_eg[i]
                                spdn_eg_ec[i] = times_ec[i] / times_eg[i]
                                spdn_ec_ic[i] = times_ic[i] / times_ec[i]

                            colors = ["cyan", "magenta", "yellow"]
                            linestyle = "-" if phase_idx == 0 else "--"
                            label_main = "spdp " + phase_label + " "
                            labels = [label_main + "ic/eg", label_main + "ec/eg", label_main + "ic/ec"]

                            axs[graph_row,graph_col].plot(vals_x, spdn_eg_ic, color=colors[0], linestyle=linestyle, label=labels[0])
                            axs[graph_row,graph_col].plot(vals_x, spdn_eg_ec, color=colors[1], linestyle=linestyle, label=labels[1])
                            axs[graph_row,graph_col].plot(vals_x, spdn_ec_ic, color=colors[2], linestyle=linestyle, label=labels[2])
                            tps[graph_row][graph_col].add_line(mytikzplot.line(vals_x, spdn_eg_ic, colors[0], linestyle, None, labels[0]))
                            tps[graph_row][graph_col].add_line(mytikzplot.line(vals_x, spdn_eg_ec, colors[1], linestyle, None, labels[1]))
                            tps[graph_row][graph_col].add_line(mytikzplot.line(vals_x, spdn_ec_ic, colors[2], linestyle, None, labels[2]))

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
                    a.set_xlabel(division_variant)
                    a.set_ylabel("Speedup")
                    tp.xlabel = division_variant
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
