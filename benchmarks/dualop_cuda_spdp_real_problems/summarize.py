#!/usr/bin/env python3

import os
import sys
import io
from datetime import datetime



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

def get_nnodes_per_side(nelems, element_type):
    if element_type == "TRIANGLE3":
        return nelems + 1
    if element_type == "TRIANGLE6":
        return 2 * nelems + 1
    if element_type == "TETRA4":
        return nelems + 1
    if element_type == "TETRA10":
        return 2 * nelems + 1
    return 0



if len(sys.argv) <= 1:
    print("not enough arguments")
    exit(1)
run_dir_name = sys.argv[1]

basedir = "benchmarks/dualop_cuda_spdp_real_problems"
run_dir = basedir + "/runs/" + run_dir_name

summdir = basedir + "/summary"
os.makedirs(summdir, exist_ok=True)

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")
summfile = summdir + "/" + datestr + ".csv"

summstring = io.StringIO()





header_1 = ["my_id", "task_id", "problem"]
header_2 = ["domain_volume", "domain_surface"]
header_3 = ["num_domains", "partition_strategy"]
header_4 = ["exitcode", "errfile", "timeoutfile"]
header_5 = ["update_time", "apply_time", "update_time_per_subdomain", "apply_time_per_subdomain"]
summstring.write(";".join(header_1))
summstring.write(";;")
summstring.write(";".join(header_2))
summstring.write(";;")
summstring.write(";".join(header_3))
summstring.write(";;")
summstring.write(";".join(header_4))
summstring.write(";;")
summstring.write(";".join(header_5))
summstring.write("\n")



tasks_dir = run_dir + "/tasks"
for task_dir_name in os.listdir(tasks_dir):
    task_dir = tasks_dir + "/" + task_dir_name
    task_id = task_dir_name

    config_file_path = task_dir + "/config.sh"
    config_file_lines = read_file_to_string(config_file_path).split("\n")

    ecf_line = [line for line in config_file_lines if ("ESPRESO_BENCHMARK_ecf_file" in line)][0]
    problem = ecf_line.split("/")[-1].split(".")[0]

    summstring.write(";")
    summstring.write(task_id + ";")
    summstring.write(problem + ";")
    summstring.write(";")
    
    out_file = task_dir + "/log_out.txt"
    out_file_lines = read_file_to_string(out_file).split("\n")

    dom_vol_lines = [l for l in out_file_lines if (" =   Domain volume [dofs]" in l)]
    if len(dom_vol_lines) == 0:
        summstring.write(";")
    else:
        line = dom_vol_lines[0]
        ndofs = line.split("[dofs]")[1].split("<")[0].replace(" ", "")
        summstring.write(ndofs + ";")

    dom_surf_lines = [l for l in out_file_lines if (" =   Domain surface [dofs]" in l)]
    if len(dom_surf_lines) == 0:
        summstring.write(";")
    else:
        line = dom_surf_lines[0]
        ndofs = line.split("[dofs]")[1].split("<")[0].replace(" ", "")
        summstring.write(ndofs + ";")

    summstring.write(";")

    line_num_dom = [line for line in config_file_lines if ("ECF_NUM_DOMAINS" in line)][0]
    value_num_doms = line_num_dom.split("=")[1]
    summstring.write(value_num_doms + ";")
    num_domains = int(value_num_doms)

    line_partstrat = [line for line in config_file_lines if ("ECF_PARTITION_STRATEGY" in line)][0]
    value_partstrat = line_partstrat.split("=")[1]
    summstring.write(value_partstrat + ";")

    summstring.write(";")

    finished_file = task_dir + "/.finished"
    exitcode = read_file_to_string(finished_file).replace("\n", "")
    summstring.write(exitcode + ";")

    err_file = task_dir + "/stderr.txt"
    err_content = read_file_to_string(err_file).replace("\n", ";").replace("\"", "\"\"")
    summstring.write("\"" + err_content + "\"" + ";")
    
    timeout_file = task_dir + "/timeout.txt"
    timeout_content = read_file_to_string(timeout_file).replace("\n", ";").replace("\"", "\"\"")
    summstring.write("\"" + timeout_content + "\"" + ";")

    summstring.write(";")

    measured_lines_update = [line for line in out_file_lines if ("finishd 'update_mainloop'" in line)]
    measured_lines_apply = [line for line in out_file_lines if ("finishd 'dualop_explicit_applicator::apply'" in line)]

    update_time = 0
    if len(measured_lines_update) > 0:
        measured_line = measured_lines_update[-1]
        update_time = float(measured_line[100:-3].replace(" ", ""))
        summstring.write(str(update_time) + ";")
    else:
        summstring.write(";")
    
    apply_time = 0
    if len(measured_lines_apply) > 2:
        measured_lines = measured_lines_apply[2:]
        apply_times = [float(l[100:-3].replace(" ", "")) for l in measured_lines]
        apply_time = sum(apply_times) / len(apply_times)
        summstring.write(str(apply_time) + ";")
    else:
        summstring.write(";")

    if len(measured_lines_update) > 0:
        avg_update_time = update_time / num_domains
        summstring.write(str(avg_update_time) + ";")
    else:
        summstring.write(";")

    if len(measured_lines_apply) > 0:
        avg_apply_time = apply_time / num_domains
        summstring.write(str(avg_apply_time) + ";")
    else:
        summstring.write(";")

    summstring.write("\n")

write_string_to_file(summfile, summstring.getvalue())
summstring.close()



print("finished. later graph with")
print("    ./benchmarks/dualop_cuda_spdp_real_problems/graphs.py " + summfile)
