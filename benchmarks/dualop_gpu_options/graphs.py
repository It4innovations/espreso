#!/usr/bin/env python3

import os
import sys
import csv
import matplotlib.pyplot as plt
import matplotlib
from datetime import datetime





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


summarize_datestr = "20240327_233112"
summarize_stage = "update"



csv_all = read_csv_to_arrays(basedir + "/summary/" + summarize_datestr + "/" + summarize_stage + ".csv")
csv_header = csv_all[0]
csv_data = csv_all[1:]

machine_col = csv_header.index("machine")
tool_col = csv_header.index("tool")
problem_col = csv_header.index("problem")
element_type_col = csv_header.index("element_type")
factor_symmetry_fake_col = csv_header.index("factor_symmetry_fake")
uniform_clusters_domains_col = csv_header.index("uniform_clusters_domains")
concurrency_update_col = csv_header.index("concurrency_update")
concurrency_apply_col = csv_header.index("concurrency_apply")
path_if_hermitian_col = csv_header.index("path_if_hermitian")
trsm1_factor_storage_col = csv_header.index("trsm1_factor_storage")
trsm2_factor_storage_col = csv_header.index("trsm2_factor_storage")
trsm1_solve_type_col = csv_header.index("trsm1_solve_type")
trsm2_solve_type_col = csv_header.index("trsm2_solve_type")
trsm_rhs_sol_order_col = csv_header.index("trsm_rhs_sol_order")
apply_scatter_gather_where_col = csv_header.index("apply_scatter_gather_where")
n_dofs_col = csv_header.index("n_dofs_per_domain")
total_col = csv_header.index("total")

csv_data0 = csv_data

problems = ["heat_transfer_2d", "heat_transfer_3d", "linear_elasticity_2d", "linear_elasticity_3d"]
for problem in problems:
    csv_data1 = list(filter(lambda row: row[problem_col] == problem, csv_data0))

    if problem[-2:] == "2d": element_types = ["TRIANGLE3", "TRIANGLE6"]
    if problem[-2:] == "3d": element_types = ["TETRA4", "TETRA10"]
    for element_type in element_types:
        csv_data2 = list(filter(lambda row: row[element_type_col] == element_type, csv_data1))

        # machines = ["karolina", "lumi"]
        # machines = ["karolina"]
        machines = ["lumi"]
        for machine in machines:
            csv_data3 = list(filter(lambda row: row[machine_col] == machine, csv_data2))

            if machine == "karolina": tools = ["cudalegacy", "cudamodern"]
            # if machine == "karolina": tools = ["cudalegacy"]
            # if machine == "karolina": tools = ["cudamodern"]
            if machine == "lumi": tools = ["rocm"]
            for tool in tools:
                plt.figure()
                fig, axs = plt.subplots(1, 2, figsize=(2560/100.0,1440/100.0))
                csv_data4 = list(filter(lambda row: row[tool_col] == tool, csv_data3))

                factor_symmetry_fake_options = ["L", "U", "B"]
                # factor_symmetry_fake_options = ["L", "U"]
                # factor_symmetry_fake_options = ["B"]
                for factor_symmetry_fake_idx in range(0,len(factor_symmetry_fake_options)):
                    factor_symmetry_fake = factor_symmetry_fake_options[factor_symmetry_fake_idx]
                    csv_data5 = list(filter(lambda row: row[factor_symmetry_fake_col] == factor_symmetry_fake, csv_data4))

                    # uniformities = ["TRUE", "FALSE"]
                    uniformities = ["TRUE"]
                    for uniform_clusters_domains in uniformities:
                        csv_data6 = list(filter(lambda row: row[uniform_clusters_domains_col] == uniform_clusters_domains, csv_data5))

                        # concurrencies_update = ["SEQ_WAIT", "SEQ_CONTINUE", "PARALLEL"]
                        concurrencies_update = ["PARALLEL"]
                        # if summarize_stage == "apply": concurrencies_update = ["PARALLEL"]
                        for concurrency_update_idx in range(0,len(concurrencies_update)):
                            concurrency_update = concurrencies_update[concurrency_update_idx]
                            csv_data7 = list(filter(lambda row: row[concurrency_update_col] == concurrency_update, csv_data6))

                            # concurrencies_apply = ["SEQ_WAIT", "SEQ_CONTINUE", "PARALLEL"]
                            if summarize_stage == "update": concurrencies_apply = [concurrency_update]
                            for concurrency_apply in concurrencies_apply:
                                csv_data8 = list(filter(lambda row: row[concurrency_apply_col] == concurrency_apply, csv_data7))

                                paths = ["TRSM", "HERK"]
                                # paths = ["TRSM"]
                                # if summarize_stage == "apply": paths = ["HERK"]
                                for path_if_hermitian_idx in range(0,len(paths)):
                                    path_if_hermitian = paths[path_if_hermitian_idx]
                                    if (factor_symmetry_fake == "L" or factor_symmetry_fake == "U") and path_if_hermitian == "TRSM": continue
                                    if (factor_symmetry_fake == "B") and path_if_hermitian == "HERK": continue
                                    csv_data9 = list(filter(lambda row: row[path_if_hermitian_col] == path_if_hermitian, csv_data8))

                                    trsm1_factor_storage_options = ["SPARSE", "DENSE"]
                                    # trsm1_factor_storage_options = ["DENSE"]
                                    # if summarize_stage == "apply": trsm1_factor_storage_options = ["SPARSE"]
                                    # if problem[-2:] == "2d": trsm1_factor_storage_options = ["SPARSE"]
                                    # if problem[-2:] == "3d": trsm1_factor_storage_options = ["DENSE"]
                                    for trsm1_factor_storage_idx in range(0,len(trsm1_factor_storage_options)):
                                        trsm1_factor_storage = trsm1_factor_storage_options[trsm1_factor_storage_idx]
                                        csv_data10 = list(filter(lambda row: row[trsm1_factor_storage_col] == trsm1_factor_storage, csv_data9))

                                        # trsm2_factor_storage_options = ["SPARSE", "DENSE"]
                                        trsm2_factor_storage_options = ["DENSE"]
                                        # if path_if_hermitian == "TRSM": trsm2_factor_storage_options = [trsm1_factor_storage]
                                        if path_if_hermitian == "HERK": trsm2_factor_storage_options = ["SPARSE"]
                                        for trsm2_factor_storage_idx in range(0,len(trsm2_factor_storage_options)):
                                            trsm2_factor_storage = trsm2_factor_storage_options[trsm2_factor_storage_idx]
                                            csv_data11 = list(filter(lambda row: row[trsm2_factor_storage_col] == trsm2_factor_storage, csv_data10))

                                            # trsm1_solve_type_options = ["L", "LHH"]
                                            trsm1_solve_type_options = [factor_symmetry_fake]
                                            if factor_symmetry_fake == "B": trsm1_solve_type_options = ["L"]
                                            # if summarize_stage == "apply": trsm1_solve_type_options = ["LHH"]
                                            for trsm1_solve_type_idx in range(0,len(trsm1_solve_type_options)):
                                                trsm1_solve_type = trsm1_solve_type_options[trsm1_solve_type_idx]
                                                csv_data12 = list(filter(lambda row: row[trsm1_solve_type_col] == trsm1_solve_type, csv_data11))

                                                # trsm2_solve_type_options = ["U", "UHH"]
                                                trsm2_solve_type_options = ["U"]
                                                # if summarize_stage == "apply": trsm2_solve_type_options = ["U"]
                                                for trsm2_solve_type_idx in range(0,len(trsm2_solve_type_options)):
                                                    trsm2_solve_type = trsm2_solve_type_options[trsm2_solve_type_idx]
                                                    csv_data13 = list(filter(lambda row: row[trsm2_solve_type_col] == trsm2_solve_type, csv_data12))

                                                    # trsm_rhs_sol_order_options = ["ROW_MAJOR", "COL_MAJOR"]
                                                    # if summarize_stage == "apply": trsm_rhs_sol_order_options = ["ROW_MAJOR"]
                                                    trsm_rhs_sol_order_options = ["ROW_MAJOR"]
                                                    # if problem[-2:] == "2d": trsm_rhs_sol_order_options = ["COL_MAJOR"]
                                                    # if problem[-2:] == "3d": trsm_rhs_sol_order_options = ["ROW_MAJOR"]
                                                    for trsm_rhs_sol_order_idx in range(0,len(trsm_rhs_sol_order_options)):
                                                        trsm_rhs_sol_order = trsm_rhs_sol_order_options[trsm_rhs_sol_order_idx]
                                                        csv_data14 = list(filter(lambda row: row[trsm_rhs_sol_order_col] == trsm_rhs_sol_order, csv_data13))

                                                        # apply_scatter_gather_where_options = ["CPU", "GPU"]
                                                        if summarize_stage == "update": apply_scatter_gather_where_options = ["GPU"]
                                                        for apply_scatter_gather_where in apply_scatter_gather_where_options:
                                                            csv_data15 = list(filter(lambda row: row[apply_scatter_gather_where_col] == apply_scatter_gather_where, csv_data14))

                                                            csv_data16 = sorted(csv_data15, key=lambda row: int(row[n_dofs_col]))
                                                            x_vals = [int(row[n_dofs_col]) for row in csv_data16]
                                                            y_vals = [(float(row[total_col]) if len(row[total_col])>0 else float("nan")) for row in csv_data16]
                                                            if len(x_vals) > 0:
                                                                linestyle = "-"
                                                                color = None
                                                                label = None
                                                                # color = "red" if ((trsm_rhs_sol_order == "ROW_MAJOR" and trsm2_solve_type == "U") or (trsm_rhs_sol_order == "COL_MAJOR" and trsm2_solve_type == "UHH")) else "blue"
                                                                # color = "red" if path_if_hermitian == "TRSM" else "blue"
                                                                # color = "red" if concurrency_update == "SEQ_WAIT" else ("green" if concurrency_update == "SEQ_CONTINUE" else "blue")
                                                                # color = "red" if trsm_rhs_sol_order == "COL_MAJOR" else "blue"
                                                                # linestyle = "-" if path_if_hermitian == "TRSM" else "--"
                                                                color = "red" if trsm1_factor_storage == "SPARSE" else "blue"
                                                                # color = "red" if trsm2_factor_storage == "SPARSE" else "blue"
                                                                # linestyle = "-" if trsm2_factor_storage == "SPARSE" else "--"
                                                                # linestyle = "-" if trsm1_factor_storage == "SPARSE" else "--"
                                                                # color = "red" if trsm1_solve_type == "L" else "blue"
                                                                # linestyle = "-" if trsm2_solve_type == "U" else "--"
                                                                # color = "red" if trsm2_solve_type == "U" else "blue"
                                                                # linestyle = "-" if trsm1_solve_type == "L" else "--"
                                                                # linestyle = "--" if factor_symmetry_fake == "L" else "-"
                                                                # color = "red" if factor_symmetry_fake == "L" else "blue"
                                                                # linestyle = "-" if trsm_rhs_sol_order == "COL_MAJOR" else "--"
                                                                # color = "red" if (factor_symmetry_fake == "L" and trsm1_solve_type == "L") or (factor_symmetry_fake == "U" and trsm1_solve_type == "LHH") else "blue"
                                                                # color = "red" if concurrency_apply == "SEQ_WAIT" else ("green" if concurrency_apply == "SEQ_CONTINUE" else "blue")
                                                                # linestyle = "-" if apply_scatter_gather_where == "CPU" else "--"
                                                                # myotherpathidx = 2 if factor_symmetry_fake == "B" else path_if_hermitian_idx
                                                                label = trsm1_factor_storage
                                                                myaxs = axs[path_if_hermitian_idx]
                                                                myaxs.loglog(x_vals, y_vals, base=2, color=color, linestyle=linestyle, label=label)

                for a in axs.flat:
                    a.legend()
                    a.set_xlabel('n_dofs')
                    a.set_ylabel('update time [ms]')
                # title = machine + "-" + tool + "-" + problem + "-" + element_type + "-" + factor_symmetry_fake + "-" + uniform_clusters_domains + "-" + concurrency_update
                title = problem + "-" + element_type + "-" + machine + "-" + tool
                # plt.xlabel('n_dofs')
                # plt.ylabel('update time [ms]')
                # plt.title(title)
                # plt.legend()
                # plt.grid(True)
                plt.savefig(graphdir + "/" + title + ".png")
                # plt.show()
                plt.close()


