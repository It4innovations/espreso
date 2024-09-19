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

def get_n_dofs_per_node(problem):
    if problem == "heat_transfer_2d":
        return 1
    elif problem == "heat_transfer_3d":
        return 1
    elif problem == "linear_elasticity_2d":
        return 2
    elif problem == "linear_elasticity_3d":
        return 3
    else:
        return -1




basedir = "benchmarks/dualop_gpu_options"
if not os.path.exists(basedir + "/summarize.py"):
    print("The script must be run from the espreso root folder")
    sys.exit(1)



# copy or put symlinks to this folder
results_dir = basedir + "/results_to_summarize"

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")
summ_dir = basedir + "/summary/" + datestr
os.makedirs(summ_dir, exist_ok=True)



runs_setupdate = 0
runs_apply = 0
for dir_name in os.listdir(results_dir):
    dir_path = results_dir + "/" + dir_name
    runs_setupdate += 2 * len(os.listdir(dir_path + "/setupdate")) + len(os.listdir(dir_path + "/setupdateapply"))
    runs_apply += len(os.listdir(dir_path + "/apply")) + len(os.listdir(dir_path + "/setupdateapply"))
runs_total = runs_setupdate + runs_apply
runs_finished = 0





outfile = summ_dir + "/set.csv"
outstring = io.StringIO()
cells_info = ["id", "str", "run_id", "machine", "tool", "dualoperator", "problem", "dofs_per_node", "element_type", "domains_x", "domains_y", "domains_z", "n_domains_total", "elements_x", "elements_y", "elements_z", "n_elements_per_domain", "avg_domain_surface", "avg_domain_volume_observed", "avg_domain_volume_calc", "factor_nnz", "concurrency_set", "concurrency_update", "concurrency_apply", "uniform_clusters_domains", "trs1_factor_storage", "trs2_factor_storage", "trs1_solve_type", "trs2_solve_type", "trsm_rhs_sol_order", "path_if_hermitian", "f_sharing_if_hermitian", "apply_scatter_gather_where", "transpose_where"]
cells_timers = ["total", "gpuinit", "gpuset", "vecresize", "sizecalc", "gpucreate", "libinit", "mainloop_outer", "mainloop_inner", "Kreg_combine", "solver_commit", "fact_symbolic", "descriptors", "buffersize", "alloc", "alloc_host", "alloc_device", "setpointers", "Bperm", "get_factors", "extract", "trans_cpu", "copyin", "kernels_preprocess", "trans_gpu", "trsm1", "trsm2", "gemm", "trsv1", "trsv2", "spmv1", "spmv2", "applystuff", "poolalloc", "wait"]
outstring.write(";".join(cells_info))
outstring.write(";;")
outstring.write("error;timeout")
outstring.write(";;")
outstring.write(";".join(cells_timers))
outstring.write("\n")
for dir_name in os.listdir(results_dir):
    dir_path = results_dir + "/" + dir_name
    subdirs = ["setupdate", "setupdateapply"]
    for subdir in subdirs:
        for run_id in os.listdir(dir_path + "/" + subdir):
            run_path = dir_path + "/" + subdir + "/" + run_id
            info_lines = read_file_to_string(run_path + "/info.txt").split("\n")
            problem = list(filter(lambda line: "ecf_file" in line, info_lines))[0].split(".")[1]
            output_lines = read_file_to_string(run_path + "/last/espreso." + problem + ".log").split("\n")
            n_domains_total = 1
            n_elements_per_domain = 1
            for field in cells_info:
                if field == "id" or field == "str" or field == "n_domains_total" or field == "n_elements_per_domain" or field == "avg_domain_volume_calc":
                    val = ""
                elif field == "run_id":
                    val = run_id
                elif field == "problem":
                    val = problem
                elif field == "dofs_per_node":
                    val = str(get_n_dofs_per_node(problem))
                elif field == "avg_domain_volume_observed":
                    lines = list(filter(lambda line: "Domain volume [dofs]" in line, output_lines))
                    if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                    else:              val = ""
                elif field == "avg_domain_surface":
                    lines = list(filter(lambda line: "Domain surface [dofs]" in line, output_lines))
                    if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                    else:              val = ""
                elif field == "factor_nnz":
                    lines = list(filter(lambda line: "Factor nnz" in line, output_lines))
                    if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                    else:              val = ""
                else:
                    val = list(filter(lambda line: field in line, info_lines))[0].split(" ")[1]
                outstring.write(val)
                outstring.write(";")

            outstring.write(";")
            outstring.write("\"" + read_file_to_string(run_path + "/last/espreso." + problem + ".err").split("\0")[0].replace("\n", ";").replace("\"", "\"\"") + "\";")
            outstring.write("\"" + read_file_to_string(run_path + "/timeout.txt").replace("\n", ";").replace("\"", "\"\"") + "\";")
            outstring.write(";")
        
            update_lines = list(filter(lambda line: "Set" in line and "rank" in line, output_lines))
            for field in cells_timers:
                lines = list(filter(lambda line: field in line, update_lines))
                if len(lines) > 0:
                    time = lines[0][70:90].replace(" ", "").replace("-nan", "").replace("nan", "")
                    outstring.write(time)
                outstring.write(";")
            outstring.write("\n")
            runs_finished += 1
            if runs_finished % 1000 == 0:
                print("Progress: ", runs_finished, "/", runs_total)
write_string_to_file(outfile, outstring.getvalue())
outstring.close()





outfile = summ_dir + "/update.csv"
outstring = io.StringIO()
cells_info = ["id", "str", "run_id", "machine", "tool", "dualoperator", "problem", "dofs_per_node", "element_type", "domains_x", "domains_y", "domains_z", "n_domains_total", "elements_x", "elements_y", "elements_z", "n_elements_per_domain", "avg_domain_surface", "avg_domain_volume_observed", "avg_domain_volume_calc", "factor_nnz", "concurrency_set", "concurrency_update", "concurrency_apply", "uniform_clusters_domains", "trs1_factor_storage", "trs2_factor_storage", "trs1_solve_type", "trs2_solve_type", "trsm_rhs_sol_order", "path_if_hermitian", "f_sharing_if_hermitian", "apply_scatter_gather_where", "transpose_where"]
cells_timers = ["total", "mainloop_outer", "mainloop_inner", "Kreg_combine", "solver_commit", "fact_numeric", "get_factors", "extract", "trans_cpu", "allocinpool", "setpointers", "copyin", "trans_gpu", "descr_update", "descr_update_trsm1", "descr_update_trsm2", "descr_update_trsv1", "descr_update_trsv2", "sp2dn", "kernels_compute", "trsm1", "trsm2", "gemm", "fcopy", "syrk", "freeinpool", "compute_d", "wait"]
outstring.write(";".join(cells_info))
outstring.write(";;")
outstring.write("error;timeout")
outstring.write(";;")
outstring.write(";".join(cells_timers))
outstring.write("\n")
for dir_name in os.listdir(results_dir):
    dir_path = results_dir + "/" + dir_name
    subdirs = ["setupdate", "setupdateapply"]
    for subdir in subdirs:
        for run_id in os.listdir(dir_path + "/" + subdir):
            run_path = dir_path + "/" + subdir + "/" + run_id
            info_lines = read_file_to_string(run_path + "/info.txt").split("\n")
            problem = list(filter(lambda line: "ecf_file" in line, info_lines))[0].split(".")[1]
            output_lines = read_file_to_string(run_path + "/last/espreso." + problem + ".log").split("\n")
            n_domains_total = 1
            n_elements_per_domain = 1
            for field in cells_info:
                if field == "id" or field == "str" or field == "n_domains_total" or field == "n_elements_per_domain" or field == "avg_domain_volume_calc":
                    val = ""
                elif field == "run_id":
                    val = run_id
                elif field == "problem":
                    val = problem
                elif field == "dofs_per_node":
                    val = str(get_n_dofs_per_node(problem))
                elif field == "avg_domain_volume_observed":
                    lines = list(filter(lambda line: "Domain volume [dofs]" in line, output_lines))
                    if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                    else:              val = ""
                elif field == "avg_domain_surface":
                    lines = list(filter(lambda line: "Domain surface [dofs]" in line, output_lines))
                    if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                    else:              val = ""
                elif field == "factor_nnz":
                    lines = list(filter(lambda line: "Factor nnz" in line, output_lines))
                    if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                    else:              val = ""
                else:
                    val = list(filter(lambda line: field in line, info_lines))[0].split(" ")[1]
                outstring.write(val)
                outstring.write(";")

            outstring.write(";")
            outstring.write("\"" + read_file_to_string(run_path + "/last/espreso." + problem + ".err").split("\0")[0].replace("\n", ";").replace("\"", "\"\"") + "\";")
            outstring.write("\"" + read_file_to_string(run_path + "/timeout.txt").replace("\n", ";").replace("\"", "\"\"") + "\";")
            outstring.write(";")
        
            update_lines = list(filter(lambda line: "Update" in line and "rank" in line, output_lines))
            for field in cells_timers:
                lines = list(filter(lambda line: field in line, update_lines))
                if len(lines) >= 2:
                    time = lines[1][70:90].replace(" ", "").replace("-nan", "").replace("nan", "")
                    outstring.write(time)
                outstring.write(";")
            outstring.write("\n")
            runs_finished += 1
            if runs_finished % 1000 == 0:
                print("Progress: ", runs_finished, "/", runs_total)
write_string_to_file(outfile, outstring.getvalue())
outstring.close()





outfile = summ_dir + "/apply.csv"
outstring = io.StringIO()
cells_info = ["id", "str", "run_id", "machine", "tool", "dualoperator", "problem", "dofs_per_node", "element_type", "domains_x", "domains_y", "domains_z", "n_domains_total", "elements_x", "elements_y", "elements_z", "n_elements_per_domain", "avg_domain_surface", "avg_domain_volume_observed", "avg_domain_volume_calc", "factor_nnz", "concurrency_set", "concurrency_update", "concurrency_apply", "uniform_clusters_domains", "trs1_factor_storage", "trs2_factor_storage", "trs1_solve_type", "trs2_solve_type", "trsm_rhs_sol_order", "path_if_hermitian", "f_sharing_if_hermitian", "apply_scatter_gather_where", "transpose_where"]
cells_timers = ["total", "copyin", "scatter", "mv_outer", "mv", "apply_outer", "allocinpool", "setpointers", "sp2dn", "compute", "spmv1", "trsv1", "trsv2", "spmv2", "freeinpool", "zerofill", "gather", "copyout", "wait"]
outstring.write(";".join(cells_info))
outstring.write(";;")
outstring.write("error;timeout")
outstring.write(";;")
outstring.write(";".join(cells_timers))
outstring.write("\n")
for dir_name in os.listdir(results_dir):
    dir_path = results_dir + "/" + dir_name
    subdirs = ["apply", "setupdateapply"]
    for subdir in subdirs:
        for run_id in os.listdir(dir_path + "/" + subdir):
            run_path = dir_path + "/" + subdir + "/" + run_id
            info_lines = read_file_to_string(run_path + "/info.txt").split("\n")
            problem = list(filter(lambda line: "ecf_file" in line, info_lines))[0].split(".")[1]
            output_lines = read_file_to_string(run_path + "/last/espreso." + problem + ".log").split("\n")
            n_domains_total = 1
            n_elements_per_domain = 1
            for field in cells_info:
                if field == "id" or field == "str" or field == "n_domains_total" or field == "n_elements_per_domain" or field == "avg_domain_volume_calc":
                    val = ""
                elif field == "run_id":
                    val = run_id
                elif field == "problem":
                    val = problem
                elif field == "dofs_per_node":
                    val = str(get_n_dofs_per_node(problem))
                elif field == "avg_domain_volume_observed":
                    lines = list(filter(lambda line: "Domain volume [dofs]" in line, output_lines))
                    if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                    else:              val = ""
                elif field == "avg_domain_surface":
                    lines = list(filter(lambda line: "Domain surface [dofs]" in line, output_lines))
                    if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                    else:              val = ""
                elif field == "factor_nnz":
                    lines = list(filter(lambda line: "Factor nnz" in line, output_lines))
                    if len(lines) > 0: val = lines[0].split("<")[0][-15:].replace(" ", "")
                    else:              val = ""
                else:
                    val = list(filter(lambda line: field in line, info_lines))[0].split(" ")[1]
                outstring.write(val)
                outstring.write(";")

            outstring.write(";")
            outstring.write("\"" + read_file_to_string(run_path + "/last/espreso." + problem + ".err").split("\0")[0].replace("\n", ";").replace("\"", "\"\"") + "\";")
            outstring.write("\"" + read_file_to_string(run_path + "/timeout.txt").replace("\n", ";").replace("\"", "\"\"") + "\";")
            outstring.write(";")
        
            update_lines = list(filter(lambda line: "Apply" in line and "rank" in line, output_lines))
            for field in cells_timers:
                lines = list(filter(lambda line: field in line, update_lines))
                if len(lines) >= 10:
                    times_str = [line[70:90] for line in lines[-10:]]
                    times = [float(tm_str) for tm_str in times_str]
                    avg = sum(times) / len(times)
                    outstring.write("{:.6f}".format(avg))
                outstring.write(";")
            outstring.write("\n")
            runs_finished += 1
            if runs_finished % 1000 == 0:
                print("Progress: ", runs_finished, "/", runs_total)
write_string_to_file(outfile, outstring.getvalue())
outstring.close()
