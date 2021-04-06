
import os
import sys

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.abspath(__file__))

class Iterator:

    def __init__(self, items):
        self.items = items
        self.keys = items.keys()
        self.pointers = [ 0 for i in items ]

    def next(self):
        self.pointers[-1] += 1
        for i in xrange(1, len(self.pointers)):
            if self.pointers[-i] == len(self.items[self.keys[-i]]):
                self.pointers[-i] = 0
                self.pointers[-i - 1] += 1

        if self.pointers[0] == len(self.items[self.keys[0]]):
            self.pointers[0] = 0
            return False

        return True

    def __getitem__(self, i):
        return self.get()[i]

    def values(self):
        return self.get().values()

    def get(self):
        result = {}
        for i in xrange(0, len(self.pointers)):
            result[self.keys[i]] = self.items[self.keys[i]][self.pointers[i]]
        return result

class ESPRESOTestEvaluator:

    def __init__(self, root):
        self.ROOT = root
        self.levels = []
        self.table = {}
        self.stats = [ "avg" ] # first, avg, min, max, ratio, count
        self.load()

    @staticmethod
    def iterate(function, *args):
        next = [ True for arg in args]
        iterators = [ Iterator(arg) for arg in args ]

        while reduce(lambda x, y: x or y, next):
            function(*[ it.get() for it in iterators ])
            for i in range(0, len(next)):
                next[i] = iterators[i].next()
                if next[i]:
                    break

    def set_stats(self, stats):
        self.stats = stats

    def get_values(self, table, position):
        if len(position):
            return self.get_values(table[position[0]], position[1:])
        else:
            return table

    def load(self):
        def fill_values(table, levels):
            if len(levels):
                for value in levels[0]:
                    table[value] = {}
                    fill_values(table[value], levels[1:])

        run = os.path.join(self.ROOT, "run")
        log = os.path.join(self.ROOT, "log")

        self.levels = []
        for path, dirs, files in os.walk(run):
            if any([file.endswith(".ecf") for file in files]):
                break
            else:
                self.levels.append(dirs)

        fill_values(self.table, self.levels)

        repetitions = max([int(file.split(".")[-2]) for file in os.listdir(log)]) + 1

        files = len(os.listdir(log))
        for file in os.listdir(log):
            position = file.split(".")[0:len(self.levels)]
            repetition = int(file.split(".")[-2])
            data = {}
            for line in open(os.path.join(log, file), "r"):
                if line.startswith("time:"):
                    measure = { i.split("=")[0].strip("'"): i.strip().split("=")[1].strip("'") for i in line[6:].split(", ") }
                    key = str(measure["region"])
                    if key not in data:
                        data[key] = { "avg": [], "min": [], "max": [], "imb": []}
                    data[key]["avg"].append(float(measure["avg"]))
                    data[key]["min"].append(float(measure["min"]))
                    data[key]["max"].append(float(measure["max"]))
                    if measure["imb"] == "> 100":
                        measure["imb"] = "100"
                    data[key]["imb"].append(float(measure["imb"]))

            values = self.get_values(self.table, position)
            for key in data:
                if key not in values:
                    values[key] = { "avg": [None] * repetitions, "min": [None] * repetitions, "max": [None] * repetitions, "imb": [None] * repetitions }
                values[key]["avg"][repetition] = sum(data[key]["avg"])
                values[key]["min"][repetition] = sum(data[key]["min"])
                values[key]["max"][repetition] = sum(data[key]["max"])
                values[key]["imb"][repetition] = sum(data[key]["imb"])


        def average(table, levels):
            if len(levels):
                for value in levels[0]:
                    average(table[value], levels[1:])
            else:
                for value in table:
                    first = table[value]["avg"][0]
                    rest = [v for v in table[value]["avg"][1:] if v is not None]
                    all = [v for v in table[value]["avg"] if v is not None]
                    amin = [v for v in table[value]["min"] if v is not None]
                    amax = [v for v in table[value]["max"] if v is not None]
                    imb = [v for v in table[value]["imb"] if v is not None]
                    if len(rest):
                        table[value] = {
                            "first": first,
                            "avg": (sum(rest) / len(rest), None)[len(rest) == 0],
                            "count": len(all),
                            "min": (sum(amin) / len(amin), None)[len(rest) == 0],
                            "max": (sum(amax) / len(amax), None)[len(rest) == 0],
                            "ratio": max(imb)}
                    else:
                        table[value] = { "first": first, "avg": None, "count": 0, "min" : None, "max": None, "ratio": None }
        average(self.table, self.levels)

    def evaluate(self, tablename, rows, columns, **kwargs):
        ranges = []
        for level in range(1, len(self.levels) + 1):
            ranges.append({"L{0}".format(level): kwargs["L{0}".format(level)]})
        ranges.append({"VALUES": kwargs["VALUES"]})

        stats = os.path.join(self.ROOT, "stats")

        trows = []
        tcolumns = []
        freevars = []
        isrowvar = []
        iscolumnvar = []
        isargfree = []
        for level in ranges:
            for k, v in level.iteritems():
                if k in rows:
                    trows.append({k: v})
                    isrowvar.append(True)
                    iscolumnvar.append(False)
                    isargfree.append(False)
                if k in columns:
                    tcolumns.append({k: v})
                    isrowvar.append(False)
                    iscolumnvar.append(True)
                    isargfree.append(False)
                if k not in rows and k not in columns:
                    freevars.append({k: v})
                    isrowvar.append(False)
                    iscolumnvar.append(False)
                    isargfree.append(True)

        tables = {}
        def add_files(tablename, stats):
            tables[tablename] = {}
            tables[tablename]["file"] = {}
            for stat in self.stats:
                tables[tablename]["file"][stat] = open(os.path.join(stats, tablename + "." + stat + ".csv"), "w")

        def create_table(*args):
            args = [ arg.values()[0].replace(" ", "_") for arg in args ]
            name = "{0}.{1}".format(tablename, ".".join(args))
            add_files(name.replace(":", ""), stats)

        if len(freevars):
            ESPRESOTestEvaluator.iterate(create_table, *freevars)
        else:
            add_files(tablename, stats)

        tablesrows = []
        tablescolumns = []
        def add_row(*args):
            tablesrows.append(".".join([ arg.values()[0].replace(" ", "_") for arg in args ]))
        def add_column(*args):
            tablescolumns.append(".".join([ arg.values()[0].replace(" ", "_") for arg in args ]))
        ESPRESOTestEvaluator.iterate(add_row, *trows)
        ESPRESOTestEvaluator.iterate(add_column, *tcolumns)

        for table in tables:
            tables[table]["data"] = {}
            for row in tablesrows:
                tables[table]["data"][row] = {}
                for column in tablescolumns:
                    tables[table]["data"][row][column] = 0;

        def write_value(*args):
            args = [ arg.values()[0] for arg in args ]
            values = self.get_values(self.table, args[:-1])
            value = None
            for key in values:
                if args[-1] == key.strip():
                    value = values[key]
                    break

            table = ".".join([tablename] + [arg.replace(" ", "_").replace(":", "") for i, arg in enumerate(args) if isargfree[i]])
            row = ".".join([arg.replace(" ", "_") for i, arg in enumerate(args) if isrowvar[i]])
            column = ".".join([arg.replace(" ", "_") for i, arg in enumerate(args) if iscolumnvar[i]])
            tables[table]["data"][row][column] = value

        ESPRESOTestEvaluator.iterate(write_value, *ranges)

        cwidth = len(max(tablescolumns, key=lambda x: len(x)))
        rwidth = len(max(tablesrows + ["-- {0} -- ".format(t) for t in tables], key=lambda x: len(x)))
        for stat in self.stats:
            for table in tables:
                tables[table]["file"][stat].write("{0:{width}} ; ".format("-- {0} --".format(table), width=rwidth))
                for column in sorted(tablescolumns):
                    tables[table]["file"][stat].write("{0:^{width}} ; ".format(column, width=cwidth))
                tables[table]["file"][stat].write("\n")
                for row in tablesrows:
                    printrow = row
                    if printrow == "SPACE":
                        printrow = ""
                    if printrow == "===":
                        printrow = "-" * rwidth
                    tables[table]["file"][stat].write("{0:<{width}} ; ".format(printrow, width=rwidth))
                    for column in sorted(tablescolumns):
                        if tables[table]["data"][row][column] is not None and tables[table]["data"][row][column][stat] is not None:
                            tables[table]["file"][stat].write("{0:>{width}.3f} ; ".format(tables[table]["data"][row][column][stat], width=cwidth))
                        else:
                            tables[table]["file"][stat].write("{0:>{width}} ; ".format("", width=cwidth))
                    tables[table]["file"][stat].write("\n")

