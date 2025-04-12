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
phase_argument = int(sys.argv[1])

basedir = "benchmarks/dualop_sctria_options"
runs_to_summarize_dir = basedir + "/runs_to_summarize"

datestr = datetime.now().strftime("%Y%m%d_%H%M%S")

summdir = basedir + "/summary/" + datestr
os.makedirs(summdir, exist_ok=True)

summfile = summdir + "/summary_phase" + str(phase_argument) + ".csv"

summstring = io.StringIO()





header_1 = ["my_id", "run_name", "machine", "task_id", "physics", "dimension"]
header_2 = ["element_type", "domains_x", "domains_y", "domains_z", "domains_total", "elements_x", "elements_y", "elements_z", "elements_total", "n_dofs", "n_dofs_domain_actual", "n_dofs_interface_actual", "dual_operator"]
header_3 = ["mainloop_update_split", "trsm_partition_algorithm", "trsm_partition_parameter", "herk_partition_algorithm", "herk_partition_parameter", "trsm_splitrhs_spdn_param", "trsm_splitfactor_gemm_spdn_param", "herk_strategy", "", "trsm_strategy", "trsm_splitrhs_spdn_criteria", "trsm_splitfactor_trsm_factor_spdn", "trsm_splitfactor_gemm_spdn_criteria", "trsm_splitfactor_gemm_factor_prune", "", "order_X", "trsm_splitrhs_factor_order_sp", "trsm_splitrhs_factor_order_dn", "trsm_splitfactor_trsm_factor_order", "trsm_splitfactor_gemm_factor_order_sp", "trsm_splitfactor_gemm_factor_order_dn"]
header_4 = ["exitcode", "errfile", "timeoutfile"]
header_5 = ["update_time", "update_time_per_subdomain", "assemble_time", "assemble_time_per_subdomain", "avg_time_trsmperform", "avg_time_herkperform"]
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


for run_dir_name in os.listdir(runs_to_summarize_dir):
    run_name = run_dir_name
    run_dir = runs_to_summarize_dir + "/" + run_dir_name
    phase_file = run_dir + "/phase.txt"
    phase_str = read_file_to_string(phase_file).replace("\n", "")
    phase_int = int(phase_str)
    if phase_int != phase_argument:
        continue
    tasks_dir = run_dir + "/tasks"
    machine_file = run_dir + "/machine.txt"
    machine = read_file_to_string(machine_file).replace("\n", "")
    for task_dir_name in os.listdir(tasks_dir):
        task_dir = tasks_dir + "/" + task_dir_name
        task_id = task_dir_name

        config_file_path = task_dir + "/config.sh"
        config_file_lines = read_file_to_string(config_file_path).split("\n")

        ecf_line = [line for line in config_file_lines if ("ESPRESO_BENCHMARK_ecf_file" in line)][0]
        problem = ecf_line.split(".")[1]
        physics = problem[0:-3]
        dimension = problem[-2:-1]

        summstring.write(";")
        summstring.write(run_name + ";")
        summstring.write(machine + ";")
        summstring.write(task_id + ";")
        summstring.write(physics + ";")
        summstring.write(dimension + ";")
        summstring.write(";")
        
        out_file = task_dir + "/out.txt"
        out_file_lines = read_file_to_string(out_file).split("\n")

        easily_extractable_fields_benchmark = ["element_type", "domains_x", "domains_y", "domains_z", "domains_total", "elements_x", "elements_y", "elements_z", "elements_total", "n_dofs", "n_dofs_domain_actual", "n_dofs_interface_actual", "dual_operator"]
        domains_total = 1
        elements_total = 1
        n_dofs = 1
        element_type = ""
        for field in easily_extractable_fields_benchmark:
            if field == "":
                summstring.write(";")
            elif field == "domains_total":
                summstring.write(str(domains_total) + ";")
            elif field == "elements_total":
                summstring.write(str(elements_total) + ";")
            elif field == "n_dofs":
                summstring.write(str(n_dofs) + ";")
            elif field == "n_dofs_domain_actual":
                lines = [l for l in out_file_lines if (" =   Domain volume [dofs]" in l)]
                if len(lines) == 0:
                    summstring.write(";")
                else:
                    line = lines[0]
                    ndofs = line.split("[dofs]")[1].split("<")[0].replace(" ", "")
                    summstring.write(ndofs + ";")
            elif field == "n_dofs_interface_actual":
                lines = [l for l in out_file_lines if (" =   Domain surface [dofs]" in l)]
                if len(lines) == 0:
                    summstring.write(";")
                else:
                    line = lines[0]
                    ndofs = line.split("[dofs]")[1].split("<")[0].replace(" ", "")
                    summstring.write(ndofs + ";")
            else:
                line = [line for line in config_file_lines if (("ESPRESO_BENCHMARK_" + field) in line)][0]
                value = line.split("\"")[1]
                summstring.write(value + ";")
                if field == "element_type":
                    element_type = value
                if "domains" in field:
                    domains_total *= int(value)
                if "elements" in field:
                    elements_total *= int(value)
                    if not (dimension == "2" and field == "elements_z"):
                        n_dofs *= get_nnodes_per_side(int(value), element_type)
        summstring.write(";")
        easily_extractable_fields_config = ["mainloop_update_split", "trsm_partition_algorithm", "trsm_partition_parameter", "herk_partition_algorithm", "herk_partition_parameter", "trsm_splitrhs_spdn_param", "trsm_splitfactor_gemm_spdn_param", "herk_strategy", "", "trsm_strategy", "trsm_splitrhs_spdn_criteria", "trsm_splitfactor_trsm_factor_spdn", "trsm_splitfactor_gemm_spdn_criteria", "trsm_splitfactor_gemm_factor_prune", "", "order_X", "trsm_splitrhs_factor_order_sp", "trsm_splitrhs_factor_order_dn", "trsm_splitfactor_trsm_factor_order", "trsm_splitfactor_gemm_factor_order_sp", "trsm_splitfactor_gemm_factor_order_dn"]
        for field in easily_extractable_fields_config:
            if field == "":
                summstring.write(";")
            else:
                line = [line for line in config_file_lines if (("ESPRESO_DUALOPSCTRIA_CONFIG_" + field) in line)][0]
                value = line.split("\"")[1]
                summstring.write(value + ";")
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

        measured_lines_wholeupdate = [line for line in out_file_lines if ("finishd 'update_mainloop'" in line)]
        measured_lines_assembleonly = [line for line in out_file_lines if ("finishd 'update_mainloop_sepatare_assemble'" in line)]

        if len(measured_lines_wholeupdate) > 0:
            measured_line = measured_lines_wholeupdate[-1]
            measured_time = measured_line[100:-3].replace(" ", "")
            summstring.write(measured_time + ";")
            time_per_subdomain = float(measured_time) / domains_total
            summstring.write(str(time_per_subdomain) + ";")
        else:
            summstring.write(";")
            summstring.write(";")

        if len(measured_lines_assembleonly) > 0:
            measured_line = measured_lines_assembleonly[-1]
            measured_time = measured_line[100:-3].replace(" ", "")
            summstring.write(measured_time + ";")
            time_per_subdomain = float(measured_time) / domains_total
            summstring.write(str(time_per_subdomain) + ";")
        else:
            summstring.write(";")
            summstring.write(";")

        times_trsmperform_cpu = [float(line[100:-3]) for line in out_file_lines if ("finishd 'trsm_csx_dny_tri::perform'" in line)]
        times_herkperform_cpu = [float(line[100:-3]) for line in out_file_lines if ("finishd 'herk_dnx_dny_tri::perform'" in line)]
        times_trsmperform_gpu = [float(line[100:-3]) for line in out_file_lines if ("finishd 'trsm_hcsx_ddny_tri::perform_submit'" in line)]
        times_herkperform_gpu = [float(line[100:-3]) for line in out_file_lines if ("finishd 'herk_ddnx_ddny_tri::perform_submit'" in line)]
        times_trsmperform_cholmod = [float(line[-12:]) for line in out_file_lines if ("tmp_trsm_cholmod_time" in line)]
        times_trsmperform_pardiso = [float(line[-15:-3]) for line in out_file_lines if ("tmp_trsm_pardiso_time" in line)]
        times_trsmperform = times_trsmperform_cpu + times_trsmperform_gpu + times_trsmperform_cholmod + times_trsmperform_pardiso
        times_herkperform = times_herkperform_cpu + times_herkperform_gpu
        if len(times_trsmperform) == 0:
            summstring.write(";")
        else:
            times_trsmperform = times_trsmperform[len(times_trsmperform)//2:]
            avg_time_trsmperform = sum(times_trsmperform) / len(times_trsmperform)
            summstring.write(str(avg_time_trsmperform) + ";")
        if len(times_herkperform) == 0:
            summstring.write(";")
        else:
            times_herkperform = times_herkperform[len(times_herkperform)//2:]
            avg_time_herkperform = sum(times_herkperform) / len(times_herkperform)
            summstring.write(str(avg_time_herkperform) + ";")

        summstring.write("\n")

write_string_to_file(summfile, summstring.getvalue())
summstring.close()
