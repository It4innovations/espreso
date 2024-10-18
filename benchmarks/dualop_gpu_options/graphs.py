#!/usr/bin/env python3

import os
import sys
import csv
import matplotlib.pyplot as plt
import matplotlib
from datetime import datetime
from mpi4py import MPI





def read_csv_to_arrays(file_path):
    rows = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file, delimiter=";")
        rows = [row for row in csv_reader]
        return rows





basedir = "benchmarks/dualop_gpu_options"
if not os.path.exists(basedir + "/graphs.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")
graphdir = basedir + "/graphs/" + datestr
os.makedirs(graphdir, exist_ok=True)


summarize_datestr = "20241015_110851"
summarize_stage = "update"



csv_all = read_csv_to_arrays(basedir + "/summary/" + summarize_datestr + "/" + summarize_stage + ".csv")
csv_header = csv_all[0]
csv_data = csv_all[1:]

machine_col = csv_header.index("machine")
tool_col = csv_header.index("tool")
problem_col = csv_header.index("problem")
element_type_col = csv_header.index("element_type")
uniform_clusters_domains_col = csv_header.index("uniform_clusters_domains")
concurrency_set_col = csv_header.index("concurrency_set")
concurrency_update_col = csv_header.index("concurrency_update")
concurrency_apply_col = csv_header.index("concurrency_apply")
path_if_hermitian_col = csv_header.index("path_if_hermitian")
trs1_factor_storage_col = csv_header.index("trs1_factor_storage")
trs2_factor_storage_col = csv_header.index("trs2_factor_storage")
trs1_solve_type_col = csv_header.index("trs1_solve_type")
trs2_solve_type_col = csv_header.index("trs2_solve_type")
trsm_rhs_sol_order_col = csv_header.index("trsm_rhs_sol_order")
apply_scatter_gather_where_col = csv_header.index("apply_scatter_gather_where")
dualoperator_col = csv_header.index("dualoperator")
n_dofs_col = csv_header.index("n_dofs")
total_col = csv_header.index("total")
factor_nnz_col = csv_header.index("factor_nnz")
# trsm1_col = csv_header.index("trsm1")
# trsm2_col = csv_header.index("trsm2")
# herk_col = csv_header.index("syrk")

csv_data0 = csv_data

#######################
subplots_counts = [12,16]
my_figsize_x = 8000
my_figsize_y = 4500
# subplots_counts = [4,12]
# my_figsize_x = 6000
# my_figsize_y = 2000
# subplots_counts = [2,3]
# my_figsize_x = 1920
# my_figsize_y = 1080
#######################
subplot_is_alone = (subplots_counts[0] == 1 and subplots_counts[1] == 1)

image_index = -1

problems = ["heat_transfer_2d", "heat_transfer_3d", "linear_elasticity_2d", "linear_elasticity_3d"]
for problem in problems:

    if problem[-2:] == "2d": element_types = ["TRIANGLE3", "TRIANGLE6"]
    if problem[-2:] == "3d": element_types = ["TETRA4", "TETRA10"]
    for element_type in element_types:
        image_index = image_index + 1
        if MPI.COMM_WORLD.Get_rank() != image_index:
            continue

        csv_data1 = list(filter(lambda row: row[problem_col] == problem, csv_data0))
        csv_data2 = list(filter(lambda row: row[element_type_col] == element_type, csv_data1))

        # machines = ["karolina", "lumi"]
        # machines = ["karolina"]
        machines = ["lumi"]
        for machine in machines:
            csv_data3 = list(filter(lambda row: row[machine_col] == machine, csv_data2))

            # tools = ["cudalegacy", "cudamodern"]
            # tools = ["cudalegacy"]
            # tools = ["cudamodern"]
            tools = ["rocm"]
            for tool in tools:
                csv_data4 = list(filter(lambda row: row[tool_col] == tool, csv_data3))

                # dualoperators = ["EXPLICIT_GPU", "IMPLICIT_GPU"]
                dualoperators = ["EXPLICIT_GPU"]
                # dualoperators = ["IMPLICIT_GPU"]
                for dualoperator in dualoperators:
                    csv_data5 = list(filter(lambda row: row[dualoperator_col] == dualoperator, csv_data4))
                    plt.figure()
                    fig, axs = plt.subplots(subplots_counts[0], subplots_counts[1], figsize=(my_figsize_x/100.0, my_figsize_y/100.0))

                    best_x_vals = ""
                    best_y_vals = ""

                    concurrencies = ["SEQ_WAIT", "SEQ_CONTINUE", "PARALLEL"]
                    # concurrencies = ["PARALLEL"]
                    # concurrencies = ["SEQ_WAIT"]
                    for concurrency_idx in range(0,len(concurrencies)):
                        concurrency = concurrencies[concurrency_idx]
                        concurrency_set = concurrency
                        concurrency_update = concurrency
                        concurrency_apply = concurrency
                        if dualoperator == "EXPLICIT_GPU" and summarize_stage == "apply": concurrency_set = "PARALLEL"
                        if dualoperator == "EXPLICIT_GPU" and summarize_stage == "apply": concurrency_update = "PARALLEL"
                        csv_data6 = list(filter(lambda row: row[concurrency_set_col]    == concurrency_set,    csv_data5))
                        csv_data7 = list(filter(lambda row: row[concurrency_update_col] == concurrency_update, csv_data6))
                        csv_data8 = list(filter(lambda row: row[concurrency_apply_col]  == concurrency_apply,  csv_data7))

                        paths = ["TRSM", "HERK"]
                        for path_idx in range(0,len(paths)):
                            path = paths[path_idx]
                            # if path != "TRSM": continue
                            # if path != "HERK": continue
                            csv_data9 = list(filter(lambda row: row[path_if_hermitian_col] == path, csv_data8))

                            trs1_factor_storage_options = ["SPARSE", "DENSE"]
                            for trs1_factor_storage_idx in range(0,len(trs1_factor_storage_options)):
                                trs1_factor_storage = trs1_factor_storage_options[trs1_factor_storage_idx]
                                csv_data10 = list(filter(lambda row: row[trs1_factor_storage_col] == trs1_factor_storage, csv_data9))
                                # if trs1_factor_storage != "DENSE": continue

                                trs2_factor_storage_options = ["SPARSE", "DENSE"]
                                # trs2_factor_storage_options = ["SPARSE"]
                                for trs2_factor_storage_idx in range(0,len(trs2_factor_storage_options)):
                                    trs2_factor_storage = trs2_factor_storage_options[trs2_factor_storage_idx]
                                    csv_data11 = list(filter(lambda row: row[trs2_factor_storage_col] == trs2_factor_storage, csv_data10))
                                    # if path == "HERK": csv_data11 = list(filter(lambda row: row[trs2_factor_storage_col] == "SPARSE", csv_data10))
                                    # if trs2_factor_storage != "DENSE": continue
                                    # if trs1_factor_storage != trs2_factor_storage: continue

                                    trs1_solve_type_options = ["L", "LHH"]
                                    # trs1_solve_type_options = ["L"]
                                    for trs1_solve_type_idx in range(0,len(trs1_solve_type_options)):
                                        trs1_solve_type = trs1_solve_type_options[trs1_solve_type_idx]
                                        csv_data12 = list(filter(lambda row: row[trs1_solve_type_col] == trs1_solve_type, csv_data11))

                                        trs2_solve_type_options = ["U", "UHH"]
                                        # trs2_solve_type_options = ["U"]
                                        # if trs1_solve_type == "L": trs2_solve_type_options = ["UHH"]
                                        # if trs1_solve_type == "LHH": trs2_solve_type_options = ["U"]
                                        for trs2_solve_type_idx in range(0,len(trs2_solve_type_options)):
                                            trs2_solve_type = trs2_solve_type_options[trs2_solve_type_idx]
                                            csv_data13 = list(filter(lambda row: row[trs2_solve_type_col] == trs2_solve_type, csv_data12))
                                            # if path == "HERK": csv_data13 = list(filter(lambda row: row[trs2_solve_type_col] == "U", csv_data12))
                                            # if trs1_solve_type == "L" and trs2_solve_type != "UHH": continue
                                            # if trs1_solve_type == "LHH" and trs2_solve_type != "U": continue

                                            trsm_rhs_sol_order_options = ["ROW_MAJOR", "COL_MAJOR"]
                                            # trsm_rhs_sol_order_options = ["ROW_MAJOR"]
                                            for trsm_rhs_sol_order_idx in range(0,len(trsm_rhs_sol_order_options)):
                                                trsm_rhs_sol_order = trsm_rhs_sol_order_options[trsm_rhs_sol_order_idx]
                                                csv_data14 = list(filter(lambda row: row[trsm_rhs_sol_order_col] == trsm_rhs_sol_order, csv_data13))
                                                # if path == "TRSM" and trsm_rhs_sol_order == "ROW_MAJOR" and not (trs1_solve_type == "L" and trs2_solve_type == "U"): continue
                                                # if path == "TRSM" and trsm_rhs_sol_order == "COL_MAJOR" and not (trs1_solve_type == "LHH" and trs2_solve_type == "UHH"): continue
                                                # if path == "HERK" and trsm_rhs_sol_order == "ROW_MAJOR" and not (trs1_solve_type == "L"): continue
                                                # if path == "HERK" and trsm_rhs_sol_order == "COL_MAJOR" and not (trs1_solve_type == "LHH"): continue
                                                # if trsm_rhs_sol_order == "ROW_MAJOR" and trs2_solve_type != "U": continue
                                                # if trsm_rhs_sol_order == "COL_MAJOR" and trs2_solve_type != "UHH": continue
                                                # if trsm_rhs_sol_order == "ROW_MAJOR" and trs1_solve_type != "L": continue
                                                # if trsm_rhs_sol_order == "COL_MAJOR" and trs1_solve_type != "LHH": continue

                                                apply_scatter_gather_where_options = ["CPU", "GPU"]
                                                # apply_scatter_gather_where_options = ["GPU"]
                                                for apply_scatter_gather_where_idx in range(0,len(apply_scatter_gather_where_options)):
                                                    apply_scatter_gather_where = apply_scatter_gather_where_options[apply_scatter_gather_where_idx]
                                                    csv_data15 = list(filter(lambda row: row[apply_scatter_gather_where_col] == apply_scatter_gather_where, csv_data14))

                                                    csv_data16 = sorted(csv_data15, key=lambda row: int(row[n_dofs_col]))
                                                    x_vals = [int(row[n_dofs_col]) for row in csv_data16]
                                                    # x_vals = [(float("nan") if row[n_dofs_col] == "" or row[factor_nnz_col] == "" else float(row[factor_nnz_col]) / float(row[n_dofs_col]) / float(row[n_dofs_col])) for row in csv_data16]
                                                    y_vals = [(float(row[total_col]) if len(row[total_col])>0 else float("nan")) for row in csv_data16]
                                                    # y_vals = [(float(row[trsm1_col]) if len(row[trsm1_col])>0 else float("nan")) for row in csv_data16]
                                                    # y_vals = [(float(row[trsm2_col]) if len(row[trsm2_col])>0 else float("nan")) for row in csv_data16]
                                                    # y_vals = [(float(row[trsm1_col])+float(row[trsm2_col]) if len(row[trsm1_col])>0 and len(row[trsm2_col])>0 else float("nan")) for row in csv_data16]
                                                    # y_vals = [(float(row[herk_col]) if len(row[herk_col])>0 else float("nan")) for row in csv_data16]
                                                    # y_vals = [(float(row[trsm1_col])+float(row[herk_col]) if len(row[trsm1_col])>0 and len(row[herk_col])>0 else float("nan")) for row in csv_data16]
                                                    if len(x_vals) > 0:
                                                        # if best_x_vals == "":
                                                        #     best_x_vals = list(x_vals)
                                                        #     best_y_vals = list(best_x_vals)
                                                        # if concurrency_apply == "SEQ_CONTINUE" and trs2_solve_type == "U":
                                                        #     if problem[-2:-1] == "2":
                                                        #         if trs1_factor_storage == "SPARSE" and trs2_factor_storage == "SPARSE" and trs1_solve_type == "L":
                                                        #             best_y_vals = list(y_vals)
                                                        #     if problem[-2:-1] == "3":
                                                        #         if trs1_factor_storage == "DENSE" and trs2_factor_storage == "DENSE" and trs1_solve_type == "LHH":
                                                        #             for j in range(0,len(x_vals)):
                                                        #                 if x_vals[j] <= 16000:
                                                        #                     best_y_vals[j] = y_vals[j]
                                                        #         if trs1_factor_storage == "SPARSE" and trs2_factor_storage == "SPARSE" and trs1_solve_type == "L":
                                                        #             for j in range(0,len(x_vals)):
                                                        #                 if x_vals[j] > 16000:
                                                        #                     best_y_vals[j] = y_vals[j]
                                                        # if concurrency_apply == "SEQ_CONTINUE" and trs1_solve_type == "L" and trs2_solve_type == "UHH":
                                                        #     if problem[-2:-1] == "2":
                                                        #         if trs1_factor_storage == "SPARSE" and trs2_factor_storage == "SPARSE":
                                                        #             best_y_vals = list(y_vals)
                                                        #     if problem[-2:-1] == "3":
                                                        #         if trs1_factor_storage == "DENSE" and trs2_factor_storage == "DENSE":
                                                        #             for j in range(0,len(x_vals)):
                                                        #                 if x_vals[j] <= 8000:
                                                        #                     best_y_vals[j] = y_vals[j]
                                                        #         if trs1_factor_storage == "SPARSE" and trs2_factor_storage == "SPARSE":
                                                        #             for j in range(0,len(x_vals)):
                                                        #                 if x_vals[j] > 8000:
                                                        #                     best_y_vals[j] = y_vals[j]
                                                        # if concurrency_apply == "PARALLEL" and trs1_solve_type == "L":
                                                        #     if problem[-2:-1] == "3":
                                                        #         if trs1_factor_storage == "DENSE" and trs2_factor_storage == "DENSE" and trs2_solve_type == "UHH":
                                                        #             best_y_vals = list(y_vals)
                                                        #     if problem[-2:-1] == "2":
                                                        #         if trs1_factor_storage == "DENSE" and trs2_factor_storage == "DENSE" and trs2_solve_type == "UHH":
                                                        #             for j in range(0,len(x_vals)):
                                                        #                 if x_vals[j] <= 8000:
                                                        #                     best_y_vals[j] = y_vals[j]
                                                        #         if trs1_factor_storage == "SPARSE" and trs2_factor_storage == "DENSE" and trs2_solve_type == "UHH":
                                                        #             for j in range(0,len(x_vals)):
                                                        #                 if x_vals[j] > 8000 and x_vals[j] <= 23000:
                                                        #                     best_y_vals[j] = y_vals[j]
                                                        #         if trs1_factor_storage == "SPARSE" and trs2_factor_storage == "SPARSE" and trs2_solve_type == "U":
                                                        #             for j in range(0,len(x_vals)):
                                                        #                 if x_vals[j] > 23000:
                                                        #                     best_y_vals[j] = y_vals[j]
                                                        # if concurrency_update == "PARALLEL" and path == "HERK" and trsm_rhs_sol_order == "ROW_MAJOR" and trs1_solve_type == "L":
                                                        #     if problem[-2:-1] == "2" and trs1_factor_storage == "SPARSE":
                                                        #         best_y_vals = list(y_vals)
                                                        #     elif problem[-2:-1] == "3":
                                                        #         if trs1_factor_storage == "SPARSE":
                                                        #             for j in range(0,len(x_vals)):
                                                        #                 if x_vals[j] > 32000:
                                                        #                     best_y_vals[j] = y_vals[j]
                                                        #         elif trs1_factor_storage == "DENSE":
                                                        #             for j in range(0,len(x_vals)):
                                                        #                 if x_vals[j] <= 32000:
                                                        #                     best_y_vals[j] = y_vals[j]
                                                        # if concurrency_update == "PARALLEL" and path == "HERK" and trs1_factor_storage == "DENSE" and trs1_solve_type == "LHH":
                                                        #     if problem[-2:-1] == "2" and trsm_rhs_sol_order == "COL_MAJOR":
                                                        #         best_y_vals = list(y_vals)
                                                        #     if problem[-2:-1] == "3" and trsm_rhs_sol_order == "ROW_MAJOR":
                                                        #         best_y_vals = list(y_vals)
                                                        # if concurrency == "PARALLEL" and trsm_rhs_sol_order == "ROW_MAJOR" and path == "HERK":
                                                        #     if problem[-2:-1] == "2" and trs1_factor_storage == "SPARSE" and trs1_solve_type == "L":
                                                        #         best_y_vals = list(y_vals)
                                                        #     if problem[-2:-1] == "3" and trs1_factor_storage == "SPARSE" and trs1_solve_type == "L":
                                                        #         for j in range(0,len(x_vals)):
                                                        #             if x_vals[j] > 12000:
                                                        #                 best_y_vals[j] = y_vals[j]
                                                        #     if problem[-2:-1] == "3" and trs1_factor_storage == "DENSE" and trs1_solve_type == "LHH":
                                                        #         for j in range(0,len(x_vals)):
                                                        #             if x_vals[j] <= 12000:
                                                        #                 best_y_vals[j] = y_vals[j]

                                                        linestyle = "-"
                                                        color = "black"
                                                        label = None
                                                        title = None
                                                        # color = "red" if concurrency == "SEQ_WAIT" else ("green" if concurrency == "SEQ_CONTINUE" else "blue")
                                                        # color = "red" if path == "TRSM" else "blue"
                                                        # linestyle = "-" if path == "TRSM" else "--"
                                                        # color = "red" if trsm_rhs_sol_order == "ROW_MAJOR" else "blue"
                                                        # linestyle = "-" if trsm_rhs_sol_order == "ROW_MAJOR" else "--"
                                                        # color = "red" if trs1_factor_storage == "SPARSE" else "blue"
                                                        # linestyle = "-" if trs1_factor_storage == "SPARSE" else "--"
                                                        # color = "red" if trs1_solve_type == "L" else "blue"
                                                        # linestyle = "-" if trs1_solve_type == "L" else "--"
                                                        # color = "red" if trs2_factor_storage == "SPARSE" else "blue"
                                                        # linestyle = "-" if trs2_factor_storage == "SPARSE" else "--"
                                                        # color = "red" if trs2_solve_type == "U" else "blue"
                                                        # linestyle = "-" if trs2_solve_type == "U" else "--"
                                                        # color = "red" if apply_scatter_gather_where == "CPU" else "blue"
                                                        # linestyle = "-" if apply_scatter_gather_where == "CPU" else "--"
                                                        color = "blue"
                                                        label = "this one"
                                                        title = concurrency + "-" + path + "-" + trsm_rhs_sol_order + "-" + trs1_factor_storage + "-" + trs2_factor_storage + "-" + trs1_solve_type + "-" + trs2_solve_type + "-" + "sg" + apply_scatter_gather_where
                                                        # title = concurrency + "-" + trs1_factor_storage + "-" + trs2_factor_storage + "-" + trs1_solve_type + "-" + trs2_solve_type
                                                        row = 4 * concurrency_idx + 2 * path_idx + trsm_rhs_sol_order_idx
                                                        col = 8 * trs1_factor_storage_idx + 4 * trs2_factor_storage_idx + 2 * trs1_solve_type_idx + trs2_solve_type_idx
                                                        # row = 2 * trs1_factor_storage_idx + trs2_factor_storage_idx
                                                        # col = 4 * concurrency_idx + 2 * trs1_solve_type_idx + trs2_solve_type_idx
                                                        # row = apply_scatter_gather_where_idx
                                                        # col = concurrency_idx
                                                        myaxs = axs[row,col]
                                                        myaxs.loglog(x_vals, y_vals, base=2, color=color, linestyle=linestyle, label=label)
                                                        if title != None: myaxs.set_title(title, fontsize="medium")

                    axs_flat = []
                    xlim_min = (1 << 30)
                    xlim_max = 0
                    ylim_min = (1 << 30)
                    ylim_max = 0
                    if subplot_is_alone: axs_flat = [axs]
                    else: axs_flat = axs.flat
                    for a in axs_flat:
                        a.grid(True)
                        # a.set_xlabel('n_dofs')
                        # a.set_ylabel('update time [ms]')
                        if a.lines:
                            # a.loglog(best_x_vals, best_y_vals, base=2, color="red", linestyle="--", label="optimal")
                            a.legend(loc="upper left")
                            xlim_min = min(xlim_min, a.get_xlim()[0])
                            xlim_max = max(xlim_max, a.get_xlim()[1])
                            ylim_min = min(ylim_min, a.get_ylim()[0])
                            ylim_max = max(ylim_max, a.get_ylim()[1])
                        a.set_xscale("log", base=2)
                        a.set_yscale("log", base=2)
                    plt.setp(axs, xlim=[xlim_min,xlim_max], ylim=[ylim_min,ylim_max])
                    title = problem + "-" + element_type + "-" + machine + "-" + tool + "-" + dualoperator
                    # plt.xlabel('n_dofs')
                    # plt.ylabel('update time [ms]')
                    # plt.title(title)
                    # plt.legend()
                    # plt.grid(True)
                    fig.tight_layout()
                    plt.savefig(graphdir + "/" + title + ".png")
                    # plt.show()
                    plt.close()



