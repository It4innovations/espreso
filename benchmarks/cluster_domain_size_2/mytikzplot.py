import random
import string





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
    def __init__(self, xs, ys, color="black", linestyle="-", marker=None, label=None, namepath=None):
        self.xs = xs
        self.ys = ys
        self.color = color
        self.linestyle = linestyle
        self.marker = marker
        self.label = label
        self.namepath = namepath
        if self.marker == None: self.marker = "none"
        if self.color == None: self.color = "only marks"
        elif self.color == "magenta": self.color = "{rgb,255:red,255;green,0;blue,255}"
        elif self.color == "cyan": self.color = "{rgb,255:red,0;green,255;blue,255}"
        elif self.color == "green": self.color = "{green!60!black}"
        elif self.color == "lawngreen": self.color = "green"
        elif self.color == "darkorange": self.color = "{orange!90!black}"
        if self.label != None: self.label = self.label.replace("_", "\\_")
    def export(self, sb):
        tikz_linestyle = "draw=none"
        if self.linestyle == "-": tikz_linestyle = "solid"
        elif self.linestyle == "--": tikz_linestyle = "dashed"
        elif self.linestyle == ":": tikz_linestyle = "dotted"
        elif self.linestyle == "-.": tikz_linestyle = "dashdotted"
        elif self.linestyle != None: print("unsupported linestyle '", self.linestyle, "'")
        forget_string = ", forget plot" if self.label == None else ""
        namepath_string = (", name path=" + self.namepath) if self.namepath != None else ""
        sb.append_line("\\addplot[" + tikz_linestyle + ", color=" + self.color + ", mark=" + self.marker + forget_string + namepath_string + "] coordinates {")
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



class area:
    def __init__(self, xs, ys_bot, ys_top, color="black", alpha=0.5, linestyle=None):
        self.xs = xs
        self.ys_bot = ys_bot
        self.ys_top = ys_top
        self.color = color
        self.alpha = alpha
        self.linestyle = linestyle
    def export(self, sb):
        name = ''.join(random.choice(string.ascii_lowercase) for _ in range(10))
        line(self.xs, self.ys_bot, color=self.color, linestyle=self.linestyle, namepath=name + "bot").export(sb)
        line(self.xs, self.ys_top, color=self.color, linestyle=self.linestyle, namepath=name + "top").export(sb)
        sb.append_line("\\addplot [" + self.color + ", fill opacity=" + str(self.alpha) + "] fill between [of=" + name + "bot and " + name + "top];")
        pass



class tikzplotter:
    def __init__(self):
        self.lines = []
        self.textboxes = []
        self.areas = []
        self.xmin = 0
        self.xmax = 100
        self.ymin = 0
        self.ymax = 100
        self.logx = 0
        self.logy = 0
        self.xlabel = None
        self.ylabel = None
        self.title = None
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
    def add_area(self, a):
        self.areas.append(a)
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
            if self.logx != 0 and self.logy == 0: axestype = "semilogyaxis"
            if self.logx == 0 and self.logy != 0: axestype = "semilogxaxis"
            if self.logx != 0 and self.logy != 0: axestype = "loglogaxis"
            sb.append_line("\\begin{" + axestype + "}[")
            sb.indent_push()
            if True:
                if self.logx != 0: sb.append_line("log basis x=" + str(self.logx) + ",")
                if self.logy != 0: sb.append_line("log basis y=" + str(self.logy) + ",")
                if self.xlabel != None: sb.append_line("xlabel={" + self.xlabel + "},")
                if self.ylabel != None: sb.append_line("ylabel={" + self.ylabel + "},")
                if self.title != None: sb.append_line("title={" + self.title + "},")
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
            for a in self.areas:
                a.export(sb)
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

