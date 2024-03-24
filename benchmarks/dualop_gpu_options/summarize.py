#!/usr/bin/env python3

import os
import sys
from datetime import datetime
import io

def read_file_to_string(file_path):
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        print("could not open file " + file_path)
        sys.exit(2)

def write_string_to_file(file_path, str):
    try:
        with open(file_path, 'w') as file:
            file.write(str)
    except FileNotFoundError:
        print("could not open file " + file_path)
        sys.exit(3)




basedir = "benchmarks/dualop_gpu_options"
if not os.path.exists(basedir + "/summarize.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)



# copy or put symlinks to this folder
results_dir = basedir + "/results_to_summarize"

# datestr = datetime.now().strftime("%Y%m%d_%H%M%S")
datestr = "20240324_232401"
summ_dir = basedir + "/summary/" + datestr
os.makedirs(summ_dir, exist_ok=True)



runs_setupdate = 0
runs_apply = 0
for dir_name in os.listdir(results_dir):
    dir_path = results_dir + "/" + dir_name
    runs_setupdate += 2 * len(os.listdir(dir_path + "/setupdate"))
    runs_apply += len(os.listdir(dir_path + "/apply"))
runs_total = runs_setupdate + runs_apply
runs_finished = 0





# outfile = summ_dir + "/set.csv"
# outstring = io.StringIO()
# cells_info = ["id", "str", "machine", "tool", "problem", "element_type", "domains_x", "domains_y", "domains_z", "n_domains_total", "elements_x", "elements_y", "elements_z", "n_elements_total", "concurrency_set", "concurrency_update", "concurrency_apply", "factor_symmetry_fake", "uniform_clusters_domains", "trsm1_factor_storage", "trsm2_factor_storage", "trsm1_solve_type", "trsm2_solve_type", "trsm_rhs_sol_order", "path_if_hermitian", "f_sharing_if_hermitian", "apply_scatter_gather_where", "transpose_where"]
# cells_timers = ["total", "gpuinit", "gpuset", "vecresize", "sizecalc", "gpucreate", "mainloop_outer", "mainloop_inner", "Kreg_combine", "solver_commit", "fact_symbolic", "descriptors", "buffersize", "alloc", "alloc_host", "alloc_device", "setpointers", "Bperm", "get_factors", "extract", "trans_cpu", "copyin", "kernels_preprocess", "trans_gpu", "trs1", "trs2", "gemm", "applystuff", "poolalloc", "wait"]
# outstring.write(";".join(cells_info))
# outstring.write(";")
# outstring.write("error;timeout")
# outstring.write(";")
# outstring.write(";".join(cells_timers))
# outstring.write("\n")
# for dir_name in os.listdir(results_dir):
#     dir_path = results_dir + "/" + dir_name
#     for run_id in os.listdir(dir_path + "/setupdate"):
#         run_path = dir_path + "/setupdate/" + run_id
#         info_lines = read_file_to_string(run_path + "/info.txt").split("\n")
#         for field in cells_info:
#             if field == "id" or field == "str" or field == "n_domains_total" or field == "n_elements_total":
#                 val = ""
#             elif field == "problem":
#                 problem = list(filter(lambda line: "ecf_file" in line, info_lines))[0].split(".")[1]
#                 val = problem
#             else:
#                 val = list(filter(lambda line: field in line, info_lines))[0].split(" ")[1]
#             outstring.write(val)
#             outstring.write(";")

#         outstring.write(";")
#         outstring.write("\"" + read_file_to_string(run_path + "/last/espreso." + problem + ".err").replace("\n", ";").replace("\"", "\"\"") + "\";")
#         outstring.write("\"" + read_file_to_string(run_path + "/timeout.txt").replace("\n", ";").replace("\"", "\"\"") + "\";")
#         outstring.write(";")
    
#         update_lines = list(filter(lambda line: "Set" in line and "rank" in line, read_file_to_string(run_path + "/last/espreso." + problem + ".log").split("\n")))
#         for field in cells_timers:
#             lines = list(filter(lambda line: field in line, update_lines))
#             if len(lines) > 0:
#                 time = lines[0][70:90].replace(" ", "")
#                 outstring.write(time)
#             outstring.write(";")
#         outstring.write("\n")
#         runs_finished += 1
#         if runs_finished % 1000 == 0:
#             print("Progress: ", runs_finished, "/", runs_total)
# write_string_to_file(outfile, outstring.getvalue())
# outstring.close()





# outfile = summ_dir + "/update.csv"
# outstring = io.StringIO()
# cells_info = ["id", "str", "machine", "tool", "problem", "element_type", "domains_x", "domains_y", "domains_z", "n_domains_total", "elements_x", "elements_y", "elements_z", "n_elements_total", "concurrency_set", "concurrency_update", "concurrency_apply", "factor_symmetry_fake", "uniform_clusters_domains", "trsm1_factor_storage", "trsm2_factor_storage", "trsm1_solve_type", "trsm2_solve_type", "trsm_rhs_sol_order", "path_if_hermitian", "f_sharing_if_hermitian", "apply_scatter_gather_where", "transpose_where"]
# cells_timers = ["total", "mainloop_outer", "mainloop_inner", "Kreg_combine", "solver_commit", "fact_numeric", "get_factors", "extract", "trans_cpu", "allocinpool", "setpointers", "copyin", "trans_gpu", "descr_update", "descr_update_trs1", "descr_update_trs2", "sp2dn", "kernels_compute", "trs1", "trs2", "gemm", "fcopy", "syrk", "freeinpool", "compute_d", "wait"]
# outstring.write(";".join(cells_info))
# outstring.write(";")
# outstring.write("error;timeout")
# outstring.write(";")
# outstring.write(";".join(cells_timers))
# outstring.write("\n")
# for dir_name in os.listdir(results_dir):
#     dir_path = results_dir + "/" + dir_name
#     for run_id in os.listdir(dir_path + "/setupdate"):
#         run_path = dir_path + "/setupdate/" + run_id
#         info_lines = read_file_to_string(run_path + "/info.txt").split("\n")
#         for field in cells_info:
#             if field == "id" or field == "str" or field == "n_domains_total" or field == "n_elements_total":
#                 val = ""
#             elif field == "problem":
#                 problem = list(filter(lambda line: "ecf_file" in line, info_lines))[0].split(".")[1]
#                 val = problem
#             else:
#                 val = list(filter(lambda line: field in line, info_lines))[0].split(" ")[1]
#             outstring.write(val)
#             outstring.write(";")

#         outstring.write(";")
#         outstring.write("\"" + read_file_to_string(run_path + "/last/espreso." + problem + ".err").replace("\n", ";").replace("\"", "\"\"") + "\";")
#         outstring.write("\"" + read_file_to_string(run_path + "/timeout.txt").replace("\n", ";").replace("\"", "\"\"") + "\";")
#         outstring.write(";")
    
#         update_lines = list(filter(lambda line: "Update" in line and "rank" in line, read_file_to_string(run_path + "/last/espreso." + problem + ".log").split("\n")))
#         for field in cells_timers:
#             lines = list(filter(lambda line: field in line, update_lines))
#             if len(lines) >= 2:
#                 time = lines[1][70:90].replace(" ", "")
#                 outstring.write(time)
#             outstring.write(";")
#         outstring.write("\n")
#         runs_finished += 1
#         if runs_finished % 1000 == 0:
#             print("Progress: ", runs_finished, "/", runs_total)
# write_string_to_file(outfile, outstring.getvalue())
# outstring.close()





outfile = summ_dir + "/apply.csv"
outstring = io.StringIO()
cells_info = ["id", "str", "machine", "tool", "problem", "element_type", "domains_x", "domains_y", "domains_z", "n_domains_total", "elements_x", "elements_y", "elements_z", "n_elements_total", "concurrency_set", "concurrency_update", "concurrency_apply", "factor_symmetry_fake", "uniform_clusters_domains", "trsm1_factor_storage", "trsm2_factor_storage", "trsm1_solve_type", "trsm2_solve_type", "trsm_rhs_sol_order", "path_if_hermitian", "f_sharing_if_hermitian", "apply_scatter_gather_where", "transpose_where"]
cells_timers = ["total", "copyin", "scatter", "mv_outer", "mv", "zerofill", "gather", "copyout", "wait"]
outstring.write(";".join(cells_info))
outstring.write(";")
outstring.write("error;timeout")
outstring.write(";")
outstring.write(";".join(cells_timers))
outstring.write("\n")
for dir_name in os.listdir(results_dir):
    dir_path = results_dir + "/" + dir_name
    for run_id in os.listdir(dir_path + "/apply"):
        run_path = dir_path + "/apply/" + run_id
        info_lines = read_file_to_string(run_path + "/info.txt").split("\n")
        for field in cells_info:
            if field == "id" or field == "str" or field == "n_domains_total" or field == "n_elements_total":
                val = ""
            elif field == "problem":
                problem = list(filter(lambda line: "ecf_file" in line, info_lines))[0].split(".")[1]
                val = problem
            else:
                val = list(filter(lambda line: field in line, info_lines))[0].split(" ")[1]
            outstring.write(val)
            outstring.write(";")

        outstring.write(";")
        outstring.write("\"" + read_file_to_string(run_path + "/last/espreso." + problem + ".err").replace("\n", ";").replace("\"", "\"\"") + "\";")
        outstring.write("\"" + read_file_to_string(run_path + "/timeout.txt").replace("\n", ";").replace("\"", "\"\"") + "\";")
        outstring.write(";")
    
        update_lines = list(filter(lambda line: "Apply" in line and "rank" in line, read_file_to_string(run_path + "/last/espreso." + problem + ".log").split("\n")))
        for field in cells_timers:
            lines = list(filter(lambda line: field in line, update_lines))
            if len(lines) >= 10:
                times_str = [line[70:90] for line in lines[-10:]]
                times = [float(tm_str) for tm_str in times_str]
                avg = sum(times) / len(times)
                outstring.write(str(avg))
            outstring.write(";")
        outstring.write("\n")
        runs_finished += 1
        if runs_finished % 1000 == 0:
            print("Progress: ", runs_finished, "/", runs_total)
write_string_to_file(outfile, outstring.getvalue())
outstring.close()
