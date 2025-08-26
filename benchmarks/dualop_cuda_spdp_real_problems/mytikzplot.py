

class my_stringbuilder:
    def __init__(self):
        self.strings = []
        self.indent = 0
        self.indent_step = 2
    def append_line(self, s):
        if self.indent > 0: self.strings.append(" " * self.indent * self.indent_step)
        self.strings.append(s)
        self.strings.append("\n")
    def indent_push(self):
        self.indent += 1
    def indent_pop(self):
        self.indent -= 1
    def get(self):
        return "".join(self.strings)





class line:
    def __init__(self, xs, ys, color, linestyle, marker, label):
        self.xs = xs
        self.ys = ys
        self.color = color
        self.linestyle = linestyle
        self.marker = marker
        self.label = label
        if self.marker == None: self.marker = "none"
        if self.color == None: self.color = "only marks"
        elif self.color == "magenta": self.color = "{rgb,255:red,255;green,0;blue,255}"
        elif self.color == "cyan": self.color = "{rgb,255:red,0;green,255;blue,255}"
        elif self.color == "green": self.color = "{green!60!black}"
        elif self.color == "lawngreen": self.color = "green"
        elif self.color == "darkorange": self.color = "{orange!90!black}"
        if self.label != None: self.label = self.label.replace("_", "\\_")
    def export(self, sb):
        if self.linestyle == "-": tikz_linestyle = "solid"
        elif self.linestyle == "--": tikz_linestyle = "dashed"
        elif self.linestyle == ":": tikz_linestyle = "dotted"
        elif self.linestyle == "-.": tikz_linestyle = "dashdotted"
        else: print("unsupported linestyle '" + self.linestyle + "'")
        forget_string = ", forget plot" if self.label == None else ""
        sb.append_line("\\addplot[" + tikz_linestyle + ", color=" + self.color + ", mark=" + self.marker + forget_string + "] coordinates {")
        sb.indent_push()
        for (x,y) in zip(self.xs,self.ys):
            sb.append_line("( {:.3e} , {:.3e} )".format(x,y))
        sb.indent_pop()
        sb.append_line("};")
            



class textbox:
    def __init__(self, text, x, y):
        self.text = text
        self.x = x
        self.y = y
    def export(self, sb):
        sb.append_line("\\node[anchor=west, inner sep=1pt, fill=white, fill opacity=0.7, text opacity=1, font=\\tiny] at (axis cs: {:.3e},{:.3e}) {{{:s}}};".format(self.x, self.y, self.text))



class tikzplotter:
    def __init__(self):
        self.lines = []
        self.textboxes = []
        self.xmin = 0
        self.xmax = 100
        self.ymin = 0
        self.ymax = 100
        self.logx = 0
        self.logy = 0
        self.xlabel = None
        self.ylabel = None
        # self.legendpos = "north west"
    @staticmethod
    def create_table_of_plotters(ny, nx):
        arr = []
        for y in range(ny):
            arr.append([])
            for x in range(nx):
                arr[y].append(tikzplotter())
        return arr
    def add_line(self, l):
        self.lines.append(l)
    def add_text(self, t):
        self.textboxes.append(t)
    # def set_legendpos(self, legendpos):
    #     self.legendpos = legendpos
    #     self.legendpos = self.legendpos.replace("upper", "north")
    #     self.legendpos = self.legendpos.replace("lower", "south")
    #     self.legendpos = self.legendpos.replace("left", "west")
    #     self.legendpos = self.legendpos.replace("right", "east")
    def set_bounds(self, xmin, xmax, ymin, ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
    def save(self, filename):
        sb = my_stringbuilder()
        sb.append_line("\\begin{tikzpicture}")
        sb.indent_push()
        if True:
            every_axis_append_styles = []
            if self.logx == 0 and self.logy == 0: axestype = "axis"
            if self.logx != 0 and self.logy == 0: axestype = "semilogxaxis"
            if self.logx == 0 and self.logy != 0: axestype = "semilogyaxis"
            if self.logx != 0 and self.logy != 0: axestype = "loglogaxis"
            sb.append_line("\\begin{" + axestype + "}[")
            sb.indent_push()
            if True:
                if self.logx != 0: sb.append_line("log basis x=" + str(self.logx) + ",")
                if self.logy != 0: sb.append_line("log basis y=" + str(self.logy) + ",")
                if self.xlabel != None: sb.append_line("xlabel={" + self.xlabel + "},")
                if self.ylabel != None: sb.append_line("ylabel={" + self.ylabel + "},")
                sb.append_line("xmajorgrids=true,")
                sb.append_line("ymajorgrids=true,")
                sb.append_line("every axis plot/.append style={line width=1.2pt},")
                sb.append_line("grid style=dotted,")
                sb.append_line("xmin={:.3e}".format(self.xmin) + ",")
                sb.append_line("xmax={:.3e}".format(self.xmax) + ",")
                sb.append_line("ymin={:.3e}".format(self.ymin) + ",")
                sb.append_line("ymax={:.3e}".format(self.ymax) + ",")
                sb.append_line("label style={font=\\footnotesize},")
                sb.append_line("tick label style={font=\\footnotesize},")
                # sb.append_line("legend pos = " + self.legendpos + ",")
                sb.append_line("legend style={font=\\footnotesize, draw=none, at={(0.5,-0.1)}, anchor=north},")
                sb.append_line("legend cell align={left}")
            sb.indent_pop()
            sb.append_line("]")
            sb.indent_push()
            for l in self.lines:
                l.export(sb)
            for t in self.textboxes:
                t.export(sb)
            legend_lines = [l for l in self.lines if l.label != None]
            for ll in legend_lines:
                sb.append_line("\\addlegendentry{" + ll.label + "}")
            sb.indent_pop()
            sb.append_line("\\end{" + axestype + "}")
        sb.indent_pop()
        sb.append_line("\\end{tikzpicture}")
        with open(filename, "w") as file:
            file.write(sb.get())

