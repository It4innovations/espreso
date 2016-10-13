
from utils import *
import unittest
import os
import glob

class ESPRESOTests(unittest.TestCase):

    espreso = Espreso()
    cube = "tests/examples/linearElasticity/cube"

    def regular_cube(self, procs, config, args):
        config["INPUT"] = "GENERATOR"
        config["PATH"] = "regular_fixed_bottom.txt"
        for file in glob.glob(self.espreso.path + "/" + self.cube + "/result*"):
            os.remove(file)
        self.espreso.run(procs, self.cube, config, args)

        vtk = 0
        vtu = 0
        vtm = 0
        for file in glob.glob(self.espreso.path + "/" + self.cube + "/result*"):
            if config["OUTPUT_FORMAT"] == "VTK_LEGACY_FORMAT":
                if file.split(".")[1] == "vtk":
                    vtk += 1
                else:
                    EspresoError("Unexpected file: {0}".format(file))
            if config["OUTPUT_FORMAT"] == "VTK_BINARY_FORMAT":
                if file.split(".")[1] == "vtu":
                    vtu += 1
                else:
                    EspresoError("Unexpected file: {0}".format(file))
            if config["OUTPUT_FORMAT"] == "VTK_MULTIBLOCK_FORMAT":
                if file.split(".")[1] == "vtu":
                    vtu += 1
                elif file.split(".")[1] == "vtm":
                    vtm += 1
                else:
                    EspresoError("Unexpected file: {0}".format(file))
            if config["OUTPUT_FORMAT"] == "ENSIGHT_FORMAT":
                print file

        if config["OUTPUT_FORMAT"] == "VTK_LEGACY_FORMAT":
            if vtk < procs:
                EspresoError("Only {0} *.vtk files for {1} procs.".format(vtk, procs))
        if config["OUTPUT_FORMAT"] == "VTK_BINARY_FORMAT":
            if vtu < procs:
                EspresoError("Only {0} *.vtu files for {1} procs.".format(vtu, procs))
        if config["OUTPUT_FORMAT"] == "VTK_MULTIBLOCK_FORMAT":
            if vtu < procs:
                EspresoError("Only {0} *.vtu files for {1} procs.".format(vtu, procs))
            if vtm != 1:
                EspresoError("Missing *.vtm file")
        if config["OUTPUT_FORMAT"] == "ENSIGHT_FORMAT":
            pass


if __name__ == '__main__':

    def parameters(config, example):
        procs = reduce(lambda x, y: x * y, example["CLUSTERS"])
        args = [example["ETYPE"]] + example["CLUSTERS"] + example["ARGS"]
        name = "_".join(str(x) for x in args + config.values())
        config["SAVE_RESULTS"] = 1
        return name, procs, args

    def regular_cube(config, example):
        name, procs, args = parameters(config, example)
        TestCaseCreator.create_test(ESPRESOTests, ESPRESOTests.regular_cube, "OUTPUT" + name, procs, config, args)

    OUTPUT_FORMATS = [ "VTK_LEGACY_FORMAT", "VTK_BINARY_FORMAT", "VTK_MULTIBLOCK_FORMAT" ] #, "ENSIGHT_FORMAT" ]
    ETYPES = [ "HEXA8", "TETRA4", "PRISMA6", "PYRAMID5", "HEXA20", "TETRA10", "PRISMA15", "PYRAMID13" ]

    # Test output format
    TestCaseCreator.iterate(
        regular_cube,
        {
            "OUTPUT_FORMAT": OUTPUT_FORMATS,
            "OUTPUT_COMPRESSION": [ 0, 1 ],
            "OUTPUT_DECIMATION": [ 0, 0.5, 0.9 ]
        },
        {
            "ETYPE": ETYPES,
            "CLUSTERS": [ [2, 2, 2] ],
            "ARGS": [ [ 2, 1, 1, 1, 2, 2] ]
        }
    )

    unittest.main()
