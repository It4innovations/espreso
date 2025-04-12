#!/usr/bin/env python3


import os
import sys
import csv
import math
import matplotlib.pyplot as plt
import matplotlib
import mytikzplot



def read_csv_to_arrays(file_path):
    rows = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file, delimiter=";")
        rows = [row for row in csv_reader]
        return rows



csv_content = read_csv_to_arrays("summ2/data.csv")
csv_header = csv_content[0]
csv_data = csv_content[1:]



col_dim = csv_header.index("dimension")
col_dualop = csv_header.index("dualop")
col_ndofs = csv_header.index("n_dofs")
col_time_update = csv_header.index("time_update")
col_time_apply = csv_header.index("time_apply")



dual_operators = [
    "suitesparse_expl",
    "pardiso_expl",
    "sctriacpu_expl",
    "sctriagpu_expl",
    "sctriagpuorig_expl",
    "suitesparse_impl",
    "pardiso_impl",
    "hybrid_expl",
    # "hypothetical",
]





csv_data_0 = csv_data
for dim in ["2", "3"]:
    csv_data_1 = [row for row in csv_data_0 if (row[col_dim] == dim)]

    imgname = "update-" + dim + "D"
    imgpath = "summ2/graphs/" + imgname + ".png"
    tikzpath = "summ2/graphs/" + imgname + ".tex"
    plt.figure()
    fig, axs = plt.subplots(1, 1, figsize=(700/100.0, 500/100.0))

    tp = mytikzplot.tikzplotter()

    for dualop in dual_operators:
        csv_data_2 = [row for row in csv_data_1 if (row[col_dualop] == dualop)]

        csv_data_3 = sorted(csv_data_2, key=lambda row: int(row[col_ndofs]))
        vals_x_str = [row[col_ndofs] for row in csv_data_3]
        vals_y_str = [row[col_time_update] for row in csv_data_3]
        vals_x = [float(x) for x in vals_x_str]
        vals_y = [(float(y) if y != "" else float("nan")) for y in vals_y_str]

        if "impl" in dualop: linestyle = "--"
        else: linestyle = "-"

        if "pardiso" in dualop: color = "cyan"
        if "suitesparse" in dualop: color = "magenta"
        if dualop == "sctriacpu_expl": color = "blue"
        if dualop == "sctriagpu_expl": color = "lawngreen"
        if dualop == "sctriagpuorig_expl": color = "green"
        if "hybrid" in dualop:
            color = "black"
            linestyle = ":"
        if dualop == "hypothetical":
            color = "black"
            linestyle = "-"

        axs.loglog(vals_x, vals_y, base=2, color=color, linestyle=linestyle, linewidth=2, label=dualop)
        tp.add_line(mytikzplot.line(vals_x, vals_y, color, linestyle, None, dualop))

    xlim_min = axs.get_xlim()[0]
    xlim_max = axs.get_xlim()[1]
    ylim_min = axs.get_ylim()[0]
    ylim_max = axs.get_ylim()[1]
    axs.grid(True)
    axs.legend()
    axs.set_xscale("log", base=2)
    axs.set_yscale("log", base=10)
    fig.tight_layout()
    plt.savefig(imgpath)
    plt.close()
    tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
    tp.logx = 2
    tp.logy = 10
    tp.xlabel = "Number of DOFs per subdomain"
    tp.ylabel = "SC assembly time per subdomain [ms]"
    tp.save(tikzpath)





csv_data_0 = csv_data
for dim in ["2", "3"]:
    csv_data_1 = [row for row in csv_data_0 if (row[col_dim] == dim)]

    imgname = "apply-" + dim + "D"
    imgpath = "summ2/graphs/" + imgname + ".png"
    tikzpath = "summ2/graphs/" + imgname + ".tex"
    plt.figure()
    fig, axs = plt.subplots(1, 1, figsize=(700/100.0, 500/100.0))

    tp = mytikzplot.tikzplotter()

    for dualop in dual_operators:
        csv_data_2 = [row for row in csv_data_1 if (row[col_dualop] == dualop)]

        csv_data_3 = sorted(csv_data_2, key=lambda row: int(row[col_ndofs]))
        vals_x_str = [row[col_ndofs] for row in csv_data_3]
        vals_y_str = [row[col_time_apply] for row in csv_data_3]
        vals_x = [float(x) for x in vals_x_str]
        vals_y = [(float(y) if y != "" else float("nan")) for y in vals_y_str]

        if "impl" in dualop: linestyle = "--"
        else: linestyle = "-"

        if "pardiso" in dualop: color = "cyan"
        if "suitesparse" in dualop: color = "magenta"
        if dualop == "sctriacpu_expl": color = "blue"
        if dualop == "sctriagpu_expl": color = "lawngreen"
        if dualop == "sctriagpuorig_expl": color = "green"
        if "hybrid" in dualop:
            color = "black"
            linestyle = ":"
        if dualop == "hypothetical":
            color = "black"
            linestyle = "-"

        axs.loglog(vals_x, vals_y, base=2, color=color, linestyle=linestyle, linewidth=2, label=dualop)
        tp.add_line(mytikzplot.line(vals_x, vals_y, color, linestyle, None, dualop))

    xlim_min = axs.get_xlim()[0]
    xlim_max = axs.get_xlim()[1]
    ylim_min = axs.get_ylim()[0]
    ylim_max = axs.get_ylim()[1]
    axs.grid(True)
    axs.legend()
    axs.set_xscale("log", base=2)
    axs.set_yscale("log", base=10)
    fig.tight_layout()
    plt.savefig(imgpath)
    plt.close()
    tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
    tp.logx = 2
    tp.logy = 10
    tp.xlabel = "Number of DOFs per subdomain"
    tp.ylabel = "Application time per subdomain [ms]"
    tp.save(tikzpath)





niters_start = 1
niters_end = 10000

csv_data_0 = csv_data
for dim in ["2", "3"]:
    csv_data_1 = [row for row in csv_data_0 if (row[col_dim] == dim)]

    imgname = "amort-" + dim + "D"
    imgpath = "summ2/graphs/" + imgname + ".png"
    tikzpath = "summ2/graphs/" + imgname + ".tex"
    plt.figure()
    fig, axs = plt.subplots(1, 1, figsize=(700/100.0, 500/100.0))

    tp = mytikzplot.tikzplotter()

    ndofs_list = sorted(list(set([int(row[col_ndofs]) for row in csv_data_1])))

    already_plotted_dualops = []

    for ndofs in ndofs_list:
        csv_data_2 = [row for row in csv_data_1 if (int(row[col_ndofs]) == ndofs)]

        curr_dualop_idx = -1
        curr_time_update = 1000000000.0
        curr_time_apply = 1000000000.0
        curr_time = 1000000000.0
        curr_niters = niters_start

        for dualop_idx in range(len(dual_operators)):
            dualop = dual_operators[dualop_idx]
            csv_data_3 = [row for row in csv_data_2 if (row[col_dualop] == dualop)]
            if len(csv_data_3) != 1:
                print("error sdfg513s4g55")
                sys.exit(1)
            row = csv_data_3[0]
            time_update_str = row[col_time_update]
            time_apply_str = row[col_time_apply]
            if len(time_update_str) == 0 or len(time_apply_str) == 0: continue
            time_update = float(time_update_str)
            time_apply = float(time_apply_str)
            time = time_update + curr_niters * time_apply
            if time < curr_time:
                curr_time = time
                curr_time_update = time_update
                curr_time_apply = time_apply
                curr_dualop_idx = dualop_idx

        while True:
            next_dualop_idx = -1
            next_niters = 1000000000.0
            next_time_update = 1000000000.0
            next_time_apply = 1000000000.0
            for dualop_idx in range(len(dual_operators)):
                dualop = dual_operators[dualop_idx]
                if dualop_idx == curr_dualop_idx: continue
                csv_data_3 = [row for row in csv_data_2 if (row[col_dualop] == dualop)]
                time_update_str = csv_data_3[0][col_time_update]
                time_apply_str = csv_data_3[0][col_time_apply]
                if len(time_update_str) == 0 or len(time_apply_str) == 0: continue
                time_update = float(time_update_str)
                time_apply = float(time_apply_str)
                if time_apply >= curr_time_apply:
                    continue
                niters = (time_update - curr_time_update) / (curr_time_apply - time_apply)
                if niters > curr_niters and niters < next_niters:
                    next_niters = niters
                    next_dualop_idx = dualop_idx
                    next_time_update = time_update
                    next_time_apply = time_apply
            if next_niters >= niters_end:
                next_niters = niters_end
            next_time = curr_time_update + next_niters * curr_time_apply

            dualop = dual_operators[curr_dualop_idx]
            if "impl" in dualop: linestyle = "--"
            else: linestyle = "-"

            if "pardiso" in dualop: color = "cyan"
            if "suitesparse" in dualop: color = "magenta"
            if dualop == "sctriacpu_expl": color = "blue"
            if dualop == "sctriagpu_expl": color = "lawngreen"
            if dualop == "sctriagpuorig_expl": color = "green"
            if "hybrid" in dualop:
                color = "black"
                linestyle = ":"
            if dualop == "hypothetical":
                color = "black"
                linestyle = "-"

            linewidth = 2
            label = dualop
            if label in already_plotted_dualops: label = None
            else: already_plotted_dualops.append(label)

            step = 1.1
            nintervals = int((math.log(next_niters) - math.log(curr_niters)) / math.log(step)) + 1
            npoints = nintervals + 1
            xs = [math.exp(x / nintervals * (math.log(next_niters) - math.log(curr_niters)) + math.log(curr_niters)) for x in range(npoints)]
            ys = [x * curr_time_apply + curr_time_update for x in xs]

            axs.loglog(xs, ys, color=color, linestyle=linestyle, linewidth=linewidth, label=label)
            tp.add_line(mytikzplot.line(xs, ys, color, linestyle, None, label))

            if next_niters >= niters_end:
                x = xs[-1] * 1.1
                y = ys[-1]
                text = str(ndofs) + " DOFs"
                axs.text(x, y, text, fontsize="medium", verticalalignment="center", bbox=dict(edgecolor="none", facecolor="white", alpha=0.7, pad=0))
                tp.add_text(mytikzplot.textbox(text, x, y))
                prev_ndofs_label = ys[-1]

            if next_niters >= niters_end:
                break
            curr_niters = next_niters
            curr_time_update = next_time_update
            curr_time_apply = next_time_apply
            curr_time = next_time
            curr_dualop_idx = next_dualop_idx

    xlim_min = axs.get_xlim()[0]
    xlim_max = axs.get_xlim()[1]
    ylim_min = axs.get_ylim()[0]
    ylim_max = axs.get_ylim()[1]
    xlim_max *= 5
    axs.grid(True)
    axs.legend(loc="upper left")
    axs.set_xscale("log", base=10)
    axs.set_yscale("log", base=10)
    plt.setp(axs, xlim=[xlim_min,xlim_max])
    fig.tight_layout()
    plt.savefig(imgpath)
    plt.close()
    tp.set_bounds(xlim_min, xlim_max, ylim_min, ylim_max)
    tp.logx = 2
    tp.logy = 10
    tp.xlabel = "Number of DOFs per subdomain"
    tp.ylabel = "Step time per subdomain [ms]"
    tp.save(tikzpath)
