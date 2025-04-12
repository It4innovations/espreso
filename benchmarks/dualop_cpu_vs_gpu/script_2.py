#!/usr/bin/env python3

import os
import io
import sys



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



outstring = io.StringIO()
header = ["machine", "dualop", "physics", "dimension", "element_type", "n_dofs", "", "time_set", "time_update", "time_apply"]
outstring.write(";".join(header))
outstring.write("\n")

topdir = "espreso_results"

for run_dir_name in os.listdir(topdir):
    run_dir_path = topdir + "/" + run_dir_name
    if not os.path.isdir(run_dir_path):
        continue
    machine = run_dir_name.split("_")[0]
    dualop = run_dir_name[14:-16]

    for task_dir_name in os.listdir(run_dir_path):
        task_dir_path = run_dir_path + "/" + task_dir_name
        infofile = task_dir_path + "/info.txt"
        outfile = task_dir_path + "/stdout.txt"
        
        info_lines = read_file_to_string(infofile).split("\n")
        out_lines = read_file_to_string(outfile).split("\n")

        problem = [row for row in info_lines if ("ecf_file" in row)][0].split(".")[-2]
        physics = problem[0:-3]
        dimension = problem[-2:-1]
        element_type = [row for row in info_lines if ("element_type" in row)][0].split(" ")[-1]

        domains_x = int([row for row in info_lines if ("domains_x" in row)][0].split(" ")[-1])
        domains_y = int([row for row in info_lines if ("domains_y" in row)][0].split(" ")[-1])
        domains_z = int([row for row in info_lines if ("domains_z" in row)][0].split(" ")[-1])
        domains_total = domains_x * domains_y * domains_z

        elements_x = int([row for row in info_lines if ("elements_x" in row)][0].split(" ")[-1])
        elements_y = int([row for row in info_lines if ("elements_y" in row)][0].split(" ")[-1])
        elements_z = int([row for row in info_lines if ("elements_z" in row)][0].split(" ")[-1])

        n_dofs = get_nnodes_per_side(elements_x, element_type) * get_nnodes_per_side(elements_y, element_type) * get_nnodes_per_side(elements_z, element_type)

        outstring.write(machine + ";")
        outstring.write(dualop + ";")
        outstring.write(physics + ";")
        outstring.write(dimension + ";")
        outstring.write(element_type + ";")
        outstring.write(str(n_dofs) + ";")

        outstring.write(";")

        times_set = [float(row.split("TMP_SET_TIME")[-1].replace(" ", "")) for row in out_lines if ("TMP_SET_TIME" in row)]
        times_update = [float(row.split("TMP_UPDATE_TIME")[-1].replace(" ", "")) for row in out_lines if ("TMP_UPDATE_TIME" in row)]
        times_apply = [float(row.split("TMP_APPLY_TIME")[-1].replace(" ", "")) for row in out_lines if ("TMP_APPLY_TIME" in row)]
        if len(times_apply) > 2:
            times_apply = times_apply[2:]

        time_set = times_set[0] / domains_total if len(times_set) > 0 else ""
        time_update = times_update[-1] / domains_total if len(times_update) > 0 else ""
        time_apply = sum(times_apply) / len(times_apply) / domains_total if len(times_apply) > 0 else ""

        outstring.write(str(time_set) + ";")
        outstring.write(str(time_update) + ";")
        outstring.write(str(time_apply))

        outstring.write("\n")

write_string_to_file("summ2/data.csv", outstring.getvalue())
outstring.close()
