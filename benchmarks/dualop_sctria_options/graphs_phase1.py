#!/usr/bin/env python3

import os
import sys
import csv
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib
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
summary_dir_timestamp = sys.argv[1]

basedir = "benchmarks/dualop_sctria_options"

summary_file = basedir + "/summary/" + summary_dir_timestamp + "/summary_phase1.csv"

if not os.path.isfile(summary_file):
    print("summary file/directory does not exit")
    exit(2)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")

graphs_dir = basedir + "/graphs/" + datestr
os.makedirs(graphs_dir, exist_ok=True)




csv_content = read_csv_to_arrays(summary_file)
csv_header = csv_content[0]
csv_data = csv_content[1:]

col_dualop = csv_header.index("dual_operator")
col_physics = csv_header.index("physics")
col_dimension = csv_header.index("dimension")
col_element_type = csv_header.index("element_type")
col_trsm_strategy = csv_header.index("trsm_strategy")
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
col_time_per_subdomain = csv_header.index("time_per_subdomain")



dualoperator_list_espreso = ["EXPLICIT_SCTRIA", "EXPLICIT_SCTRIA_GPU"]
dualoperator_list_cpugpu = ["CPU", "GPU"]
physics_list = ["heat_transfer", "linear_elasticity"]
dimensions_list = ["2", "3"]
element_types_2d_list = ["TRIANGLE3", "TRIANGLE6"]
element_types_3d_list = ["TETRA4", "TETRA10"]

trsm_strategy_list = ["R", "F", "_"]
trsm_splitrhs_spdn_list = ["S", "D", "_"]
trsm_splitfactor_trsm_spdn_list = ["S", "D", "_"]
trsm_splitfactor_gemm_spdn_list = ["S", "D", "_"]
trsm_splitfactor_gemm_prune_list = ["N", "R", "_"]

order_X_list = ["R", "C", "_"]
trsm_splitrhs_factor_order_sp_list = ["R", "C", "_"]
trsm_splitrhs_factor_order_dn_list = ["R", "C", "_"]
trsm_splitfactor_trsm_factor_order_list = ["R", "C", "_"]
trsm_splitfactor_gemm_factor_order_sp_list = ["R", "C", "_"]
trsm_splitfactor_gemm_factor_order_dn_list = ["R", "C", "_"]

subplots_counts = [8,8]
my_figsize_x = 4000
my_figsize_y = 3000

image_index = -1

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
                for trsm_strategy in trsm_strategy_list:
                    csv_data_05 = [row for row in csv_data_04 if (row[col_trsm_strategy] == trsm_strategy)]
                    for trsm_splitrhs_spdn in trsm_splitrhs_spdn_list:
                        csv_data_06 = [row for row in csv_data_05 if (row[col_trsm_splitrhs_spdn] == trsm_splitrhs_spdn)]
                        for trsm_splitfactor_trsm_spdn in trsm_splitfactor_trsm_spdn_list:
                            csv_data_07 = [row for row in csv_data_06 if (row[col_trsm_splitfactor_trsm_spdn] == trsm_splitfactor_trsm_spdn)]
                            for trsm_splitfactor_gemm_spdn in trsm_splitfactor_gemm_spdn_list:
                                csv_data_08 = [row for row in csv_data_07 if (row[col_trsm_splitfactor_gemm_spdn] == trsm_splitfactor_gemm_spdn)]
                                for trsm_splitfactor_gemm_prune in trsm_splitfactor_gemm_prune_list:
                                    csv_data_09 = [row for row in csv_data_08 if(row[col_trsm_splitfactor_gemm_prune] == trsm_splitfactor_gemm_prune)]
                                    if len(csv_data_09) == 0:
                                        continue
                                    image_index += 1
                                    if MPI.COMM_WORLD.Get_rank() != image_index:
                                        continue

                                    imgname = dualop_cpugpu + "-" + physics + "-" + dimension + "D-" + element_type + "-" + trsm_strategy + trsm_splitrhs_spdn + trsm_splitfactor_trsm_spdn + trsm_splitfactor_gemm_spdn + trsm_splitfactor_gemm_prune
                                    imgpath = graphs_dir + "/" + imgname + ".png"
                                    plt.figure()
                                    fig, axs = plt.subplots(subplots_counts[0], subplots_counts[1], figsize=(my_figsize_x/100.0, my_figsize_y/100.0))

                                    for order_X_idx in range(len(order_X_list)):
                                        order_X = order_X_list[order_X_idx]
                                        csv_data_10 = [row for row in csv_data_09 if (row[col_order_X] == order_X)]
                                        for trsm_splitrhs_factor_order_sp_idx in range(len(trsm_splitrhs_factor_order_sp_list)):
                                            trsm_splitrhs_factor_order_sp = trsm_splitrhs_factor_order_sp_list[trsm_splitrhs_factor_order_sp_idx]
                                            csv_data_11 = [row for row in csv_data_10 if (row[col_trsm_splitrhs_factor_order_sp] == trsm_splitrhs_factor_order_sp)]
                                            for trsm_splitrhs_factor_order_dn_idx in range(len(trsm_splitrhs_factor_order_dn_list)):
                                                trsm_splitrhs_factor_order_dn = trsm_splitrhs_factor_order_dn_list[trsm_splitrhs_factor_order_dn_idx]
                                                csv_data_12 = [row for row in csv_data_11 if (row[col_trsm_splitrhs_factor_order_dn] == trsm_splitrhs_factor_order_dn)]
                                                for trsm_splitfactor_trsm_factor_order_idx in range(len(trsm_splitfactor_trsm_factor_order_list)):
                                                    trsm_splitfactor_trsm_factor_order = trsm_splitfactor_trsm_factor_order_list[trsm_splitfactor_trsm_factor_order_idx]
                                                    csv_data_13 = [row for row in csv_data_12 if (row[col_trsm_splitfactor_trsm_factor_order] == trsm_splitfactor_trsm_factor_order)]
                                                    for trsm_splitfactor_gemm_factor_order_sp_idx in range(len(trsm_splitfactor_gemm_factor_order_sp_list)):
                                                        trsm_splitfactor_gemm_factor_order_sp = trsm_splitfactor_gemm_factor_order_sp_list[trsm_splitfactor_gemm_factor_order_sp_idx]
                                                        csv_data_14 = [row for row in csv_data_13 if (row[col_trsm_splitfactor_gemm_factor_order_sp] == trsm_splitfactor_gemm_factor_order_sp)]
                                                        for trsm_splitfactor_gemm_factor_order_dn_idx in range(len(trsm_splitfactor_gemm_factor_order_dn_list)):
                                                            trsm_splitfactor_gemm_factor_order_dn = trsm_splitfactor_gemm_factor_order_dn_list[trsm_splitfactor_gemm_factor_order_dn_idx]
                                                            csv_data_15 = [row for row in csv_data_14 if (row[col_trsm_splitfactor_gemm_factor_order_dn] == trsm_splitfactor_gemm_factor_order_dn)]
                                                            if(len(csv_data_15) == 0):
                                                                continue

                                                            csv_data_16 = sorted(csv_data_15, key=lambda row: int(row[col_n_dofs]))
                                                            vals_x_str = [row[col_n_dofs] for row in csv_data_16]
                                                            vals_y_str = [row[col_time_per_subdomain] for row in csv_data_16]
                                                            vals_x = [float(x) for x in vals_x_str]
                                                            vals_y = [(float(y) if y != "" else float("nan")) for y in vals_y_str]

                                                            use_order_X_idx = order_X_idx % 2
                                                            use_trsm_splitrhs_factor_order_sp_idx = trsm_splitrhs_factor_order_sp_idx % 2
                                                            use_trsm_splitrhs_factor_order_dn_idx = trsm_splitrhs_factor_order_dn_idx % 2
                                                            use_trsm_splitfactor_trsm_factor_order_idx = trsm_splitfactor_trsm_factor_order_idx % 2
                                                            use_trsm_splitfactor_gemm_factor_order_sp_idx = trsm_splitfactor_gemm_factor_order_sp_idx % 2
                                                            use_trsm_splitfactor_gemm_factor_order_dn_idx = trsm_splitfactor_gemm_factor_order_dn_idx % 2

                                                            graph_row = 4 * use_trsm_splitrhs_factor_order_sp_idx + 2 * use_trsm_splitfactor_gemm_factor_order_sp_idx + use_trsm_splitfactor_trsm_factor_order_idx
                                                            graph_col = 4 * use_trsm_splitrhs_factor_order_dn_idx + 2 * use_trsm_splitfactor_gemm_factor_order_dn_idx + use_order_X_idx

                                                            color = "black"
                                                            linestyle = ":"
                                                            # if AAA == "R":
                                                            #     color = "red"
                                                            #     linestyle = "-"
                                                            # if AAA == "C":
                                                            #     color = "blue"
                                                            #     linestyle = "--"
                                                            label = "this one"
                                                            # label = AAA
                                                            label = label.replace("_", "-")
                                                            title = "RTs" + trsm_splitrhs_factor_order_sp + "-FGs" + trsm_splitfactor_gemm_factor_order_sp + "-FT" + trsm_splitfactor_trsm_factor_order + "---RTd" + trsm_splitrhs_factor_order_dn + "-FGd" + trsm_splitfactor_gemm_factor_order_dn + "-X" + order_X
                                                            myaxs = axs[graph_row,graph_col]
                                                            myaxs.loglog(vals_x, vals_y, base=2, color=color, linestyle=linestyle, label=label)
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













