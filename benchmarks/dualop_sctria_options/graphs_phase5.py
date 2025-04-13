#!/usr/bin/env python3

import os
import sys
import csv
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib
from mpi4py import MPI
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
summary_dir_timestamp = sys.argv[1]

basedir = "benchmarks/dualop_sctria_options"

summary_file = basedir + "/summary/" + summary_dir_timestamp + "/summary_phase5.csv"

if not os.path.isfile(summary_file):
    print("summary file/directory does not exit")
    exit(2)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")

graphs_dir = basedir + "/graphs/" + datestr
os.makedirs(graphs_dir, exist_ok=True)
os.makedirs(graphs_dir + "/tikz", exist_ok=True)




csv_content = read_csv_to_arrays(summary_file)
csv_header = csv_content[0]
csv_data = csv_content[1:]

col_dualop = csv_header.index("dual_operator")
col_physics = csv_header.index("physics")
col_dimension = csv_header.index("dimension")
col_element_type = csv_header.index("element_type")
col_trsm_strategy = csv_header.index("trsm_strategy")
col_trsm_partition_parameter = csv_header.index("trsm_partition_parameter")
col_trsm_splitrhs_spdn = csv_header.index("trsm_splitrhs_spdn_criteria")
col_trsm_splitfactor_trsm_spdn = csv_header.index("trsm_splitfactor_trsm_factor_spdn")
col_trsm_splitfactor_gemm_spdn = csv_header.index("trsm_splitfactor_gemm_spdn_criteria")
col_trsm_splitfactor_gemm_prune = csv_header.index("trsm_splitfactor_gemm_factor_prune")
col_order_X = csv_header.index("order_X")
col_trsm_splitrhs_factor_order_sp = csv_header.index("trsm_splitrhs_factor_order_sp")
col_trsm_splitrhs_factor_order_dn = csv_header.index("trsm_splitrhs_factor_order_dn")
col_trsm_splitfactor_trsm_factor_order = csv_header.index("trsm_splitfactor_trsm_factor_order")
col_trsm_splitfactor_gemm_factor_order_sp = csv_header.index("trsm_splitfactor_gemm_factor_order_sp")
col_trsm_splitfactor_gemm_factor_order_dn = csv_header.index("trsm_splitfactor_gemm_factor_order_dn")
col_n_dofs = csv_header.index("n_dofs")
col_time_per_subdomain = csv_header.index("assemble_time_per_subdomain")
col_herk_strategy = csv_header.index("herk_strategy")
col_herk_partition_parameter = csv_header.index("herk_partition_parameter")
col_trsm_partition_algorithm = csv_header.index("trsm_partition_algorithm")
col_herk_partition_algorithm = csv_header.index("herk_partition_algorithm")



dualoperator_list_espreso = ["EXPLICIT_SCTRIA", "EXPLICIT_SCTRIA_GPU"]
dualoperator_list_cpugpu = ["CPU", "GPU"]
physics_list = ["heat_transfer", "linear_elasticity"]
dimensions_list = ["2", "3"]
element_types_2d_list = ["TRIANGLE3", "TRIANGLE6"]
element_types_3d_list = ["TETRA4", "TETRA10"]

trsm_strategy_list = ["R", "F"]
trsm_partition_parameter_list_2d_interface = ["-500", "-50", "10", "100"]
trsm_partition_parameter_list_3d_interface = ["-1000", "-100", "10", "100"]
trsm_partition_parameter_list_2d_domain = ["-20000", "-2000", "10", "100"]
trsm_partition_parameter_list_3d_domain = ["-2000", "-200", "10", "100"]
trsm_splitrhs_spdn_list = ["D", "S", "_"]
trsm_splitfactor_trsm_spdn_list = ["D", "S", "_"]
trsm_splitfactor_gemm_spdn_list = ["D", "S", "_"]
trsm_splitfactor_gemm_prune_list = ["N", "R", "_"]
herk_strategy_list = ["T", "Q"]
trsm_partition_algorithm_list = ["U", "M"]
herk_partition_algorithm_list = ["U", "M"]

image_index = -1





subplots_counts = [2,4]
my_figsize_x = 3000
my_figsize_y = 1000

csv_data_00 = csv_data
for dualop_idx in range(len(dualoperator_list_espreso)):
    dualop = dualoperator_list_espreso[dualop_idx]
    dualop_cpugpu = dualoperator_list_cpugpu[dualop_idx]
    csv_data_01 = [row for row in csv_data_00 if (row[col_dualop] == dualop)]
    for physics in physics_list:
        csv_data_02 = [row for row in csv_data_01 if (row[col_physics] == physics)]
        for dimension in dimensions_list:
            csv_data_03 = [row for row in csv_data_02 if (row[col_dimension] == dimension)]
            if dimension == "2": element_types_list = element_types_2d_list
            if dimension == "3": element_types_list = element_types_3d_list
            for element_type in element_types_list:
                csv_data_04 = [row for row in csv_data_03 if (row[col_element_type] == element_type)]
                if(len(csv_data_04) == 0):
                    continue
                image_index += 1
                if MPI.COMM_WORLD.Get_rank() != image_index:
                    continue

                imgname = dualop_cpugpu + "-" + physics + "-" + dimension + "D-" + element_type + "-trsm"
                imgpath = graphs_dir + "/" + imgname + ".png"
                plt.figure()
                fig, axs = plt.subplots(subplots_counts[0], subplots_counts[1], figsize=(my_figsize_x/100.0, my_figsize_y/100.0))

                tps = mytikzplot.tikzplotter.create_table_of_plotters(subplots_counts[0], subplots_counts[1])

                for trsm_partition_algorithm_idx in range(len(trsm_partition_algorithm_list)):
                    trsm_partition_algorithm = trsm_partition_algorithm_list[trsm_partition_algorithm_idx]
                    csv_data_05 = [row for row in csv_data_04 if (row[col_trsm_partition_algorithm] == trsm_partition_algorithm)]
                    herk_partition_algorithm = "_"
                    csv_data_06 = [row for row in csv_data_05 if (row[col_herk_partition_algorithm] == herk_partition_algorithm)]
                    for trsm_strategy_idx in range(len(trsm_strategy_list)):
                        trsm_strategy = trsm_strategy_list[trsm_strategy_idx]
                        csv_data_07 = [row for row in csv_data_06 if (row[col_trsm_strategy] == trsm_strategy)]
                        herk_strategy = "_"
                        csv_data_08 = [row for row in csv_data_07 if (row[col_herk_strategy] == herk_strategy)]
                        for trsm_splitfactor_gemm_prune_idx in range(len(trsm_splitfactor_gemm_prune_list)):
                            trsm_splitfactor_gemm_prune = trsm_splitfactor_gemm_prune_list[trsm_splitfactor_gemm_prune_idx]
                            csv_data_09 = [row for row in csv_data_08 if (row[col_trsm_splitfactor_gemm_prune] == trsm_splitfactor_gemm_prune)]
                            if(len(csv_data_09) == 0):
                                continue

                            csv_data_10 = sorted(csv_data_09, key=lambda row: int(row[col_n_dofs]))
                            vals_x_str = [row[col_n_dofs] for row in csv_data_10]
                            vals_y_str = [row[col_time_per_subdomain] for row in csv_data_10]
                            vals_x = [float(x) for x in vals_x_str]
                            vals_y = [(float(y) if y != "" else float("nan")) for y in vals_y_str]

                            graph_row = trsm_splitfactor_gemm_prune_idx
                            graph_col = 2 * trsm_partition_algorithm_idx + trsm_strategy_idx

                            color = "black"
                            linestyle = ":"
                            # color = "blue" if trsm_strategy == "R" else "red"
                            # linestyle = "--" if trsm_splitfactor_gemm_prune == "R" else "-"
                            # if trsm_strategy == "R":
                            #     color = "blue"
                            #     linestyle = "-"
                            # if trsm_strategy == "F":
                            #     color = "red"
                            #     linestyle = "--"
                            label = "this one"
                            # label = trsm_strategy + "-" + trsm_splitfactor_gemm_prune
                            label = label.replace("_", "-")
                            title = "partalg" + trsm_partition_algorithm + "-strategy" + trsm_strategy + "-prune" + trsm_splitfactor_gemm_prune
                            myaxs = axs[graph_row,graph_col]
                            myaxs.loglog(vals_x, vals_y, base=2, color=color, linestyle=linestyle, label=label)
                            tps[graph_row][graph_col].add_line(mytikzplot.line(vals_x, vals_y, color, linestyle, None, label))
                            # tps[graph_row][graph_col].add_line(mytikzplot.line(vals_x, vals_y, color, linestyle, None, label))
                            if title != None: myaxs.set_title(title, fontsize="medium")

                xlim_min = (1 << 30)
                xlim_max = 0
                ylim_min = (1 << 30)
                ylim_max = 0
                for a in axs.flat:
                    a.grid(True)
                    if a.lines:
                        a.legend(loc="upper left")
                        xlim_min = min(xlim_min, a.get_xlim()[0])
                        xlim_max = max(xlim_max, a.get_xlim()[1])
                        ylim_min = min(ylim_min, a.get_ylim()[0])
                        ylim_max = max(ylim_max, a.get_ylim()[1])
                    a.set_xscale("log", base=2)
                    a.set_yscale("log", base=2)
                plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
                title = imgname
                # plt.title(title)
                fig.tight_layout()
                plt.savefig(imgpath)
                plt.close()

                for trsm_partition_algorithm_idx in range(len(trsm_partition_algorithm_list)):
                    trsm_partition_algorithm = trsm_partition_algorithm_list[trsm_partition_algorithm_idx]
                    for trsm_strategy_idx in range(len(trsm_strategy_list)):
                        trsm_strategy = trsm_strategy_list[trsm_strategy_idx]
                        for trsm_splitfactor_gemm_prune_idx in range(len(trsm_splitfactor_gemm_prune_list)):
                            trsm_splitfactor_gemm_prune = trsm_splitfactor_gemm_prune_list[trsm_splitfactor_gemm_prune_idx]
                            if trsm_splitfactor_gemm_prune == "_": continue
                            graph_row = trsm_splitfactor_gemm_prune_idx
                            graph_col = 2 * trsm_partition_algorithm_idx + trsm_strategy_idx
                            tp = tps[graph_row][graph_col]
                            tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
                            tp.logx = 2
                            tp.logy = 10
                            tp.xlabel = "Number of DOFs per subdomain"
                            tp.ylabel = "SC assembly time per subdomain [ms]"
                            title = "partalg" + trsm_partition_algorithm + "-strategy" + trsm_strategy + "-prune" + trsm_splitfactor_gemm_prune
                            tikzpath = graphs_dir + "/tikz/" + imgname + "-" + title + ".tex"
                            tp.save(tikzpath)






subplots_counts = [2,2]
my_figsize_x = 1500
my_figsize_y = 1000

csv_data_00 = csv_data
for dualop_idx in range(len(dualoperator_list_espreso)):
    dualop = dualoperator_list_espreso[dualop_idx]
    dualop_cpugpu = dualoperator_list_cpugpu[dualop_idx]
    csv_data_01 = [row for row in csv_data_00 if (row[col_dualop] == dualop)]
    for physics in physics_list:
        csv_data_02 = [row for row in csv_data_01 if (row[col_physics] == physics)]
        for dimension in dimensions_list:
            csv_data_03 = [row for row in csv_data_02 if (row[col_dimension] == dimension)]
            if dimension == "2": element_types_list = element_types_2d_list
            if dimension == "3": element_types_list = element_types_3d_list
            for element_type in element_types_list:
                csv_data_04 = [row for row in csv_data_03 if (row[col_element_type] == element_type)]
                if(len(csv_data_04) == 0):
                    continue
                image_index += 1
                if MPI.COMM_WORLD.Get_rank() != image_index:
                    continue

                imgname = dualop_cpugpu + "-" + physics + "-" + dimension + "D-" + element_type + "-herk"
                imgpath = graphs_dir + "/" + imgname + ".png"
                plt.figure()
                fig, axs = plt.subplots(subplots_counts[0], subplots_counts[1], figsize=(my_figsize_x/100.0, my_figsize_y/100.0))

                tps = mytikzplot.tikzplotter.create_table_of_plotters(subplots_counts[0], subplots_counts[1])

                trsm_partition_algorithm = "_"
                csv_data_05 = [row for row in csv_data_04 if (row[col_trsm_partition_algorithm] == trsm_partition_algorithm)]
                for herk_partition_algorithm_idx in range(len(herk_partition_algorithm_list)):
                    herk_partition_algorithm = herk_partition_algorithm_list[herk_partition_algorithm_idx]
                    csv_data_06 = [row for row in csv_data_05 if (row[col_herk_partition_algorithm] == herk_partition_algorithm)]
                    trsm_strategy = "_"
                    csv_data_07 = [row for row in csv_data_06 if (row[col_trsm_strategy] == trsm_strategy)]
                    for herk_strategy_idx in range(len(herk_strategy_list)):
                        herk_strategy = herk_strategy_list[herk_strategy_idx]
                        csv_data_08 = [row for row in csv_data_07 if (row[col_herk_strategy] == herk_strategy)]
                        trsm_splitfactor_gemm_prune = "_"
                        csv_data_09 = [row for row in csv_data_08 if (row[col_trsm_splitfactor_gemm_prune] == trsm_splitfactor_gemm_prune)]
                        if(len(csv_data_09) == 0):
                            continue

                        csv_data_10 = sorted(csv_data_09, key=lambda row: int(row[col_n_dofs]))
                        vals_x_str = [row[col_n_dofs] for row in csv_data_10]
                        vals_y_str = [row[col_time_per_subdomain] for row in csv_data_10]
                        vals_x = [float(x) for x in vals_x_str]
                        vals_y = [(float(y) if y != "" else float("nan")) for y in vals_y_str]

                        graph_row = herk_partition_algorithm_idx
                        graph_col = herk_strategy_idx

                        color = "black"
                        linestyle = ":"
                        # if herk_strategy == "T":
                        #     color = "blue"
                        #     linestyle = "-"
                        # if herk_strategy == "Q":
                        #     color = "red"
                        #     linestyle = "--"
                        label = "this one"
                        # label = herk_strategy
                        label = label.replace("_", "-")
                        title = "partalg" + herk_partition_algorithm + "-strategy" + herk_strategy
                        myaxs = axs[graph_row,graph_col]
                        myaxs.loglog(vals_x, vals_y, base=2, color=color, linestyle=linestyle, label=label)
                        tps[graph_row][graph_col].add_line(mytikzplot.line(vals_x, vals_y, color, linestyle, None, label))
                        if title != None: myaxs.set_title(title, fontsize="medium")

                xlim_min = (1 << 30)
                xlim_max = 0
                ylim_min = (1 << 30)
                ylim_max = 0
                for a in axs.flat:
                    a.grid(True)
                    if a.lines:
                        a.legend(loc="upper left")
                        xlim_min = min(xlim_min, a.get_xlim()[0])
                        xlim_max = max(xlim_max, a.get_xlim()[1])
                        ylim_min = min(ylim_min, a.get_ylim()[0])
                        ylim_max = max(ylim_max, a.get_ylim()[1])
                    a.set_xscale("log", base=2)
                    a.set_yscale("log", base=2)
                plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
                title = imgname
                # plt.title(title)
                fig.tight_layout()
                plt.savefig(imgpath)
                plt.close()

                for herk_partition_algorithm_idx in range(len(herk_partition_algorithm_list)):
                    herk_partition_algorithm = herk_partition_algorithm_list[herk_partition_algorithm_idx]
                    for herk_strategy_idx in range(len(herk_strategy_list)):
                        herk_strategy = herk_strategy_list[herk_strategy_idx]
                        graph_row = herk_partition_algorithm_idx
                        graph_col = herk_strategy_idx
                        tp = tps[graph_row][graph_col]
                        tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
                        tp.logx = 2
                        tp.logy = 10
                        tp.xlabel = "Number of DOFs per subdomain"
                        tp.ylabel = "SC assembly time per subdomain [ms]"
                        title = "partalg" + herk_partition_algorithm + "-strategy" + herk_strategy
                        tikzpath = graphs_dir + "/tikz/" + imgname + "-" + title + ".tex"
                        tp.save(tikzpath)












