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

summary_file = basedir + "/summary/" + summary_dir_timestamp + "/summary_phase7.csv"

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
col_assemble_time_per_subdomain = csv_header.index("assemble_time_per_subdomain")
col_update_time_per_subdomain = csv_header.index("update_time_per_subdomain")
col_herk_strategy = csv_header.index("herk_strategy")
col_herk_partition_parameter = csv_header.index("herk_partition_parameter")
col_mainloop_update_split = csv_header.index("mainloop_update_split")
col_avg_time_trsmperform = csv_header.index("avg_time_trsmperform")
col_avg_time_herkperform = csv_header.index("avg_time_herkperform")



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
mainloop_update_split_list = ["S", "C"]

subplots_counts = [2,2]
my_figsize_x = 1500
my_figsize_y = 1000

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
                if(len(csv_data_04) == 0):
                    continue
                image_index += 1
                if MPI.COMM_WORLD.Get_rank() != image_index:
                    continue

                imgname = dualop_cpugpu + "-" + physics + "-" + dimension + "D-" + element_type
                imgpath = graphs_dir + "/" + imgname + ".png"
                plt.figure()
                fig, axs = plt.subplots(subplots_counts[0], subplots_counts[1], figsize=(my_figsize_x/100.0, my_figsize_y/100.0))

                csv_data_05 = csv_data_04
                csv_data_06 = sorted(csv_data_05, key=lambda row: int(row[col_n_dofs]))

                csv_data_baseline = [row for row in csv_data_06 if (row[col_trsm_partition_parameter] == "1")]
                csv_data_optimal = [row for row in csv_data_06 if (row[col_trsm_partition_parameter] == "0")]

                vals_x_str = [row[col_n_dofs] for row in csv_data_baseline]
                vals_x = [float(x) for x in vals_x_str]

                vals_y_trsm_baseline_str = [row[col_avg_time_trsmperform] for row in csv_data_baseline]
                vals_y_trsm_optimal_str = [row[col_avg_time_trsmperform] for row in csv_data_optimal]
                vals_y_trsm_baseline = [(float(y) if y != "" else float("nan")) for y in vals_y_trsm_baseline_str]
                vals_y_trsm_optimal = [(float(y) if y != "" else float("nan")) for y in vals_y_trsm_optimal_str]
                vals_y_speedup_trsm = [float("nan")] * len(vals_y_trsm_baseline)
                for i in range(len(vals_y_speedup_trsm)):
                    vals_y_speedup_trsm[i] = vals_y_trsm_baseline[i] / vals_y_trsm_optimal[i]

                vals_y_herk_baseline_str = [row[col_avg_time_herkperform] for row in csv_data_baseline]
                vals_y_herk_optimal_str = [row[col_avg_time_herkperform] for row in csv_data_optimal]
                vals_y_herk_baseline = [(float(y) if y != "" else float("nan")) for y in vals_y_herk_baseline_str]
                vals_y_herk_optimal = [(float(y) if y != "" else float("nan")) for y in vals_y_herk_optimal_str]
                vals_y_speedup_herk = [float("nan")] * len(vals_y_herk_baseline)
                for i in range(len(vals_y_speedup_herk)):
                    vals_y_speedup_herk[i] = vals_y_herk_baseline[i] / vals_y_herk_optimal[i]

                axs[0,0].loglog(vals_x, vals_y_trsm_baseline, color="blue", linestyle="-", label="baseline time")
                axs[0,0].loglog(vals_x, vals_y_trsm_optimal, color="red", linestyle="--", label="optimized time")
                axs[1,0].semilogx(vals_x, vals_y_speedup_trsm, color="green", linestyle="-", label="speedup")

                axs[0,1].loglog(vals_x, vals_y_herk_baseline, color="blue", linestyle="-", label="baseline time")
                axs[0,1].loglog(vals_x, vals_y_herk_optimal, color="red", linestyle="--", label="optimized time")
                axs[1,1].semilogx(vals_x, vals_y_speedup_herk, color="green", linestyle="-", label="speedup")

                axs[0,0].set_title("individual trsm kernel", fontsize="medium")
                axs[1,0].set_title("individual trsm kernel", fontsize="medium")
                axs[0,1].set_title("individual herk kernel", fontsize="medium")
                axs[1,1].set_title("individual herk kernel", fontsize="medium")


                xlim_min = (1 << 30)
                xlim_max = 0
                for graph_row in range(2):
                    ylim_min = (1 << 30)
                    ylim_max = 0
                    for graph_col in range(2):
                        a = axs[graph_row, graph_col]
                        a.grid(True)
                        if a.lines:
                            a.legend(loc="upper left")
                            xlim_min = min(xlim_min, a.get_xlim()[0])
                            xlim_max = max(xlim_max, a.get_xlim()[1])
                            ylim_min = min(ylim_min, a.get_ylim()[0])
                            ylim_max = max(ylim_max, a.get_ylim()[1])
                        a.set_xscale("log", base=2)
                        if graph_row == 0: a.set_yscale("log", base=2)
                    if graph_row == 0:
                        for graph_col in range(2):
                            a = axs[graph_row, graph_col]
                            a.set_ylim([ylim_min,ylim_max])
                plt.setp(axs, xlim=[xlim_min,xlim_max])
                fig.tight_layout()
                plt.savefig(imgpath)
                plt.close()

















