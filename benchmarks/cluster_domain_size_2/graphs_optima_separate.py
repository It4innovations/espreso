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
graphs_dir_png = graphs_dir + "/png_optima_separate"
graphs_dir_tikz = graphs_dir + "/tikz_optima_separate"
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

bound_1 = 1.2
bound_2 = 1.4
bound_3 = 2

image_index = -1

csv_data_00 = csv_data
for physics in physics_list:
    for element_type in element_types_list:
        for approach in approach_list:
            image_index = image_index + 1
            if (image_index % MPI.COMM_WORLD.Get_size()) != MPI.COMM_WORLD.Get_rank():
                continue

            csv_data_03 = [row for row in csv_data_00 if ((row[col_physics] == physics) and (row[col_element_type] == element_type) and (row[col_approach] == approach))]
            dimension = csv_data_03[0][col_dimension]

            n_dofs_cluster_list = sorted(list(set([int(row[col_n_dofs_cluster]) for row in csv_data_03])))

            graph_nrows = len(division_variant_cols) * len(part_labels)
            graph_ncols = len(feti_method_list) * len(phase_labels)

            imgname = physics + "-" + dimension + "D-" + element_type + "-" + approach
            imgpath = graphs_dir_png + "/" + imgname + ".png"
            plt.figure()
            fig, axs = plt.subplots(graph_nrows, graph_ncols, figsize=(my_figsize_x*graph_ncols/100.0, my_figsize_y*graph_nrows/100.0), sharex=True, sharey="row")

            tps = mytikzplot.tikzplotter.create_table_of_plotters(graph_nrows, graph_ncols)

            for feti_method_idx in range(len(feti_method_list)):
                feti_method = feti_method_list[feti_method_idx]
                csv_data_04 = [row for row in csv_data_03 if (row[col_feti_method] == feti_method)]

                for division_variant_idx in range(len(division_variant_cols)):
                    division_variant_col = division_variant_cols[division_variant_idx]
                    csv_data_05 = sorted(csv_data_04, key=lambda row: int(row[division_variant_col]))

                    for part_idx in range(len(part_labels)):
                        part = part_labels[part_idx]

                        for phase_idx in range(len(phase_labels)):
                            phase = phase_labels[phase_idx]
                            value_col = part_phase_cols[part_idx][phase_idx]

                            graph_row = part_idx + len(part_labels) * division_variant_idx
                            graph_col = phase_idx + len(phase_labels) * feti_method_idx

                            title = feti_method + "-" + part + "-" + phase
                            axs[graph_row,graph_col].set_title(title)
                            tps[graph_row][graph_col].title = title

                            vals_x = [int(x) for x in n_dofs_cluster_list]
                            vals_y_highest = [float("nan")] * len(vals_x)
                            vals_y_higher  = [float("nan")] * len(vals_x)
                            vals_y_high    = [float("nan")] * len(vals_x)
                            vals_y_optimum = [float("nan")] * len(vals_x)
                            vals_y_low     = [float("nan")] * len(vals_x)
                            vals_y_lower   = [float("nan")] * len(vals_x)
                            vals_y_lowest  = [float("nan")] * len(vals_x)

                            for n_dofs_cluster_idx in range(len(n_dofs_cluster_list)):
                                n_dofs_cluster = n_dofs_cluster_list[n_dofs_cluster_idx]
                                csv_data_06 = [row for row in csv_data_05 if (row[col_n_dofs_cluster] == str(n_dofs_cluster))]
                                times_str = [row[value_col] for row in csv_data_06]
                                times = [(float("nan") if len(x) == 0 else float(x)) for x in times_str]
                                division_variants_str = [row[division_variant_col] for row in csv_data_06]
                                division_variants = [(float("nan") if len(x) == 0 else float(x)) for x in division_variants_str]

                                index_optimum = -1
                                optimum_val = float("inf")
                                for i,v in enumerate(times):
                                    if v == v:
                                        if v < optimum_val:
                                            optimum_val = v
                                            index_optimum = i
                                if index_optimum == -1:
                                    continue

                                index_high = index_optimum
                                index_low = index_optimum
                                while index_high < len(times)-1 and times[index_high] < bound_1 * times[index_optimum]: index_high += 1
                                while index_low > 0 and times[index_low] < bound_1 * times[index_optimum]: index_low -= 1
                                index_higher = index_high
                                index_lower = index_low
                                while index_higher < len(times)-1 and times[index_higher] < bound_2 * times[index_optimum]: index_higher += 1
                                while index_lower > 0 and times[index_lower] < bound_2 * times[index_optimum]: index_lower -= 1
                                index_highest = index_higher
                                index_lowest = index_lower
                                while index_highest < len(times)-1 and times[index_highest] < bound_3 * times[index_optimum]: index_highest += 1
                                while index_lowest > 0 and times[index_lowest] < bound_3 * times[index_optimum]: index_lowest -= 1

                                vals_y_highest[n_dofs_cluster_idx] = division_variants[index_highest]
                                vals_y_higher [n_dofs_cluster_idx] = division_variants[index_higher]
                                vals_y_high   [n_dofs_cluster_idx] = division_variants[index_high]
                                vals_y_optimum[n_dofs_cluster_idx] = division_variants[index_optimum]
                                vals_y_low    [n_dofs_cluster_idx] = division_variants[index_low]
                                vals_y_lower  [n_dofs_cluster_idx] = division_variants[index_lower]
                                vals_y_lowest [n_dofs_cluster_idx] = division_variants[index_lowest]

                            color = "red" if feti_method == "TOTAL_FETI" else "blue"
                            linestyle = "-"
                            axs[graph_row,graph_col].fill_between(vals_x, vals_y_lowest, vals_y_highest, color=color, alpha=0.2, linewidth=0)
                            axs[graph_row,graph_col].fill_between(vals_x, vals_y_lower, vals_y_higher, color=color, alpha=0.2, linewidth=0)
                            axs[graph_row,graph_col].fill_between(vals_x, vals_y_low, vals_y_high, color=color, alpha=0.2, linewidth=0)
                            axs[graph_row,graph_col].plot(vals_x, vals_y_optimum, color=color, linestyle=linestyle, label=None)
                            tps[graph_row][graph_col].add_area(mytikzplot.area(vals_x, vals_y_lowest, vals_y_highest, color=color, alpha=0.2))
                            tps[graph_row][graph_col].add_area(mytikzplot.area(vals_x, vals_y_lower, vals_y_higher, color=color, alpha=0.2))
                            tps[graph_row][graph_col].add_area(mytikzplot.area(vals_x, vals_y_low, vals_y_high, color=color, alpha=0.2))
                            tps[graph_row][graph_col].add_line(mytikzplot.line(vals_x, vals_y_optimum, color, linestyle, None, None))

            for graph_row in range(graph_nrows):
                for graph_col in range(graph_ncols):
                    tp = tps[graph_row][graph_col]
                    a = axs[graph_row, graph_col]
                    a.grid(True)
                    a.xaxis.set_tick_params(labelbottom=True)
                    a.yaxis.set_tick_params(labelleft=True)
                    a.set_xscale("log", base=10)
                    a.set_yscale("log", base=10)
                    tp.logx = 10
                    tp.logy = 10
                    a.set_xlabel("Number of DOFs in cluster")
                    a.set_ylabel(division_variant_labels[graph_row // len(part_labels)])
                    tp.xlabel = "Number of DOFs in cluster"
                    tp.ylabel = division_variant_labels[graph_row // len(part_labels)]
                    [xlim_min, xlim_max] = a.get_xlim()
                    [ylim_min, ylim_max] = a.get_ylim()
                    tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
            fig.tight_layout()
            plt.savefig(imgpath)
            plt.close()

            for graph_row in range(graph_nrows):
                for graph_col in range(graph_ncols):
                    path = graphs_dir_tikz + "/" + imgname + "-" + tps[graph_row][graph_col].title + "-" + division_variant_labels_lite[graph_row // len(part_labels)] + ".tex"
                    tps[graph_row][graph_col].save(path)

