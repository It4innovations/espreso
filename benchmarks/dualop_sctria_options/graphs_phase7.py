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



dualoperator_list_espreso = ["EXPLICIT_SCTRIA", "EXPLICIT_SCTRIA_GPU", "EXPLICIT_SC_CHOLMOD", "EXPLICIT_PARDISO"]
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
for cpugpu_idx in range(len(dualoperator_list_cpugpu)):
    dualop_cpugpu = dualoperator_list_cpugpu[cpugpu_idx]
    csv_data_01 = [row for row in csv_data_00 if (("GPU" in row[col_dualop]) == ("GPU" in dualop_cpugpu))]
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

                csv_data_orig = [row for row in csv_data_06 if (row[col_trsm_partition_parameter] == "1")]
                csv_data_sctria = [row for row in csv_data_06 if (row[col_trsm_partition_parameter] == "0")]
                csv_data_cholmod = [row for row in csv_data_06 if (row[col_dualop] == "EXPLICIT_SC_CHOLMOD")]
                csv_data_pardiso = [row for row in csv_data_06 if (row[col_dualop] == "EXPLICIT_PARDISO")]

                vals_x_str = [row[col_n_dofs] for row in csv_data_orig]
                vals_x = [float(x) for x in vals_x_str]

                vals_y_trsm_orig_str = [row[col_avg_time_trsmperform] for row in csv_data_orig]
                vals_y_trsm_sctria_str = [row[col_avg_time_trsmperform] for row in csv_data_sctria]
                vals_y_trsm_orig = [(float(y) if y != "" else float("nan")) for y in vals_y_trsm_orig_str]
    
                vals_y_trsm_sctria = [(float(y) if y != "" else float("nan")) for y in vals_y_trsm_sctria_str]
                vals_y_speedup_sctria = [float("nan")] * len(vals_y_trsm_orig)
                for i in range(len(vals_y_speedup_sctria)):
                    vals_y_speedup_sctria[i] = vals_y_trsm_orig[i] / vals_y_trsm_sctria[i]
                
                if len(csv_data_cholmod) > 0:
                    vals_y_trsm_cholmod_str = [row[col_avg_time_trsmperform] for row in csv_data_cholmod]
                    vals_y_trsm_cholmod = [(float(y) if y != "" else float("nan")) for y in vals_y_trsm_cholmod_str]
                    vals_y_speedup_cholmod = [float("nan")] * len(vals_y_trsm_orig)
                    for i in range(len(vals_y_speedup_cholmod)):
                        vals_y_speedup_cholmod[i] = vals_y_trsm_cholmod[i] / vals_y_trsm_sctria[i]
                if len(csv_data_pardiso) > 0:
                    vals_y_trsm_pardiso_str = [row[col_avg_time_trsmperform] for row in csv_data_pardiso]
                    vals_y_trsm_pardiso = [(float(y) if y != "" else float("nan")) for y in vals_y_trsm_pardiso_str]
                    vals_y_speedup_pardiso = [float("nan")] * len(vals_y_trsm_orig)
                    for i in range(len(vals_y_speedup_pardiso)):
                        vals_y_speedup_pardiso[i] = vals_y_trsm_pardiso[i] / vals_y_trsm_sctria[i]

                vals_y_herk_orig_str = [row[col_avg_time_herkperform] for row in csv_data_orig]
                vals_y_herk_sctria_str = [row[col_avg_time_herkperform] for row in csv_data_sctria]
                vals_y_herk_orig = [(float(y) if y != "" else float("nan")) for y in vals_y_herk_orig_str]
                vals_y_herk_sctria = [(float(y) if y != "" else float("nan")) for y in vals_y_herk_sctria_str]
                vals_y_speedup_herk = [float("nan")] * len(vals_y_herk_orig)
                for i in range(len(vals_y_speedup_herk)):
                    vals_y_speedup_herk[i] = vals_y_herk_orig[i] / vals_y_herk_sctria[i]

                axs[0,0].loglog(vals_x, vals_y_trsm_orig, color="green", linestyle="-", label="orig time")
                axs[0,0].loglog(vals_x, vals_y_trsm_sctria, color="lawngreen", linestyle="-", label="sctria time")
                axs[1,0].loglog(vals_x, vals_y_speedup_sctria, color="green", linestyle="-", label="speedup over orig")
                if len(csv_data_cholmod) > 0:
                    axs[0,0].loglog(vals_x, vals_y_trsm_cholmod, color="magenta", linestyle="-", label="cholmod time")
                    axs[1,0].loglog(vals_x, vals_y_speedup_cholmod, color="magenta", linestyle="-", label="speedup over cholmod")
                if len(csv_data_pardiso) > 0:
                    axs[0,0].loglog(vals_x, vals_y_trsm_pardiso, color="cyan", linestyle="-", label="pardiso time")
                    axs[1,0].loglog(vals_x, vals_y_speedup_pardiso, color="cyan", linestyle="-", label="speedup over pardiso")

                axs[0,1].loglog(vals_x, vals_y_herk_orig, color="green", linestyle="-", label="basorigeline time")
                axs[0,1].loglog(vals_x, vals_y_herk_sctria, color="lawngreen", linestyle="-", label="sctria time")
                axs[1,1].loglog(vals_x, vals_y_speedup_herk, color="green", linestyle="-", label="speedup over orig")

                axs[0,0].set_title("individual trsm kernel, time", fontsize="medium")
                axs[1,0].set_title("individual trsm kernel, speedup of sctria over X", fontsize="medium")
                axs[0,1].set_title("individual herk kernel, time", fontsize="medium")
                axs[1,1].set_title("individual herk kernel, speedup of sctria over X", fontsize="medium")


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
                        if graph_row == 0: a.set_yscale("log", base=10)
                        if graph_row == 1: a.set_yscale("log", base=10)
                    if graph_row == 0:
                        for graph_col in range(2):
                            a = axs[graph_row, graph_col]
                            a.set_ylim([ylim_min,ylim_max])
                plt.setp(axs, xlim=[xlim_min,xlim_max])
                fig.tight_layout()
                plt.savefig(imgpath)
                plt.close()

















