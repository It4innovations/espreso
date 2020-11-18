
import sys, os, logging, subprocess, types

from waflib import Logs
from waflib.Build import BuildContext
from waflib.Configure import ConfigurationContext

class ShowConfiguration(BuildContext):
    cmd = "show"
    fun = "show"

from waflib.Build import BuildContext
class ShowEnv(BuildContext):
    cmd = "env"
    fun = "env"

def link_cxx(self, *k, **kw):
    includes = []
    libpath = []
    if "root" in kw and os.path.isdir(kw["root"]):
        includes = [ os.path.join(kw["root"], dir) for dir in os.listdir(kw["root"]) if dir.startswith("include") ]
        libpath = [ os.path.join(kw["root"], dir) for dir in os.listdir(kw["root"]) if dir.startswith("lib") ]

    general = dict(uselib_store=kw["name"].upper(), mandatory=False)
    if "use" in kw:
        general["use"] = kw["use"]

    header = dict()
    if "header_name" in kw:
        header = dict(header_name=kw["header_name"], define_name="", defines=["HAVE_" + kw["name"].upper()], includes=includes)
        header.update(general)
        header["msg"] = "Checking for {0} header ".format(kw["name"])
        if not self.check_cxx(**header):
            return False

        if "fragment" in kw:
            self.env["DEFINES_" + kw["name"].upper()] = []
            self.env["INCLUDES_" + kw["name"].upper()] = []
            test = dict(execute=True)
            test.update(header)
            inc = [ "#include <{0}>\n".format(h) for h in kw["header_name"].split() ]
            test["fragment"] = "{0}int main(int argc, char** argv) {{ {1} }}".format("".join(inc), kw["fragment"])
            test["msg"] = "Checking for {0} settings ".format(kw["name"])
            if not self.check_cxx(**test):
                return False

        if "libs" in kw:
            # remove defines since libs are not checked yet
            self.env["DEFINES_" + kw["name"].upper()] = []
            self.env["INCLUDES_" + kw["name"].upper()] = []

    libs = dict(stlib=kw["libs"], libpath=libpath, msg="Checking for {0} static library ".format(kw["name"]))
    libs.update(general)
    libs.update(header)
    if not self.options.static or not self.check_cxx(**libs):
        libs["lib"] = libs["stlib"]
        libs.pop("stlib")
        libs["msg"] = "Checking for {0} dynamic library ".format(kw["name"])
        return self.check_cxx(**libs)

    return True

from waflib.TaskGen import after_method, feature
@feature('cxx')
@after_method('apply_link')
def reorder_libs(ctx):
    def myorder(lib):
        if lib == "parmetis":
            return 1
        if lib == "metis":
            return 2
        if lib == "pardiso":
            return 3
        if lib.startswith("mkl"):
            return 4
        return 0

    ctx.env["LIB"].sort(key=myorder)

def print_available(ctx):
    def _print(msg, err, libs, color="RED"):
        libs = [lib for lib in libs if "HAVE_" + lib.upper() in ctx.env["DEFINES_" + lib.upper()]]
        ctx.start_msg(msg)
        if len(libs) == 0:
            ctx.end_msg(err, color=color)
        else:
            ctx.end_msg("[ " + ", ".join(libs) + " ]", color="BLUE")
        return bool(len(libs))

    ctx.env["HAVE_DECOMPOSERS"] = _print(
        "Available graph partitioning tools",
        "NOT FOUND [ESPRESO FUNCTIONALITY IS SIGNIFICANTLY LIMITED]",
        [ "metis", "parmetis", "scotch", "ptscotch", "kahip" ])
    ctx.env["HAVE_MATH"] = _print(
        "Available math libraries",
        "NOT FOUND [ESPRESO SOLVER CANNOT BE COMPILED]",
        [ "mkl" ])
    _print(
        "Available third party solvers",
        "NOT FOUND",
        [ "mklpdss", "hypre", "pardiso", "superlu", "wsmp" ],
        "YELLOW")
    _print(
        "Available miscellaneous libraries",
        "NOT FOUND",
        [ "async", "hdf5", "bem", "catayst" ],
        "YELLOW")

""" Recurse to third party libraries wrappers"""
def recurse(ctx):
    """ Graph partition tools """
    ctx.recurse("src/wrappers/metis")
    ctx.recurse("src/wrappers/parmetis")
    ctx.recurse("src/wrappers/scotch")
    ctx.recurse("src/wrappers/ptscotch")
    ctx.recurse("src/wrappers/kahip")

    """ Math libraries"""
    ctx.recurse("src/wrappers/mkl")

    """ Solvers """
    ctx.recurse("src/wrappers/mklpdss")
    ctx.recurse("src/wrappers/hypre")
    ctx.recurse("src/wrappers/pardiso")
    ctx.recurse("src/wrappers/superlu")
    ctx.recurse("src/wrappers/wsmp")

    """ Other """
    ctx.recurse("src/wrappers/async")
    ctx.recurse("src/wrappers/hdf5")
    ctx.recurse("src/wrappers/bem")
    ctx.recurse("src/wrappers/catalyst")

def configure(ctx):
    ctx.link_cxx = types.MethodType(link_cxx, ctx)
    def trycompiler():
        try:
            ctx.find_program(ctx.options.mpicxx, var="MPI_CXX")
        except ctx.errors.ConfigurationError:
            return False
        return True

    """ Set compilers """

    if not trycompiler():
        if ctx.options.mpicxx == "mpiicpc":
            ctx.options.mpicxx = "mpic++"
        elif ctx.options.mpicxx == "mpic++":
            ctx.options.mpicxx = "mpiicpc"
        if not trycompiler():
            ctx.fatal("Cannot found MPI compiler. Set a correct one by 'mpicxx=' parameter.")

    if ctx.options.mpicxx == "mpiicpc":
        ctx.options.cxx = "icpc"
    if ctx.options.mpicxx == "mpic++":
        ctx.options.cxx = "g++"
    ctx.load(ctx.options.cxx)
    ctx.env.CXX = ctx.env.LINK_CXX = ctx.env.MPI_CXX

    """ Set default compilers flags"""

    ctx.env.append_unique("CXXFLAGS", [ "-fopenmp" ])
    ctx.env.append_unique("LINKFLAGS", [ "-fopenmp" ])

    if ctx.options.intwidth == "32":
        ctx.env.append_unique("DEFINES", [ "esint=int", "esint_mpi=MPI_INT" ])
    if ctx.options.intwidth == "64":
        ctx.env.append_unique("DEFINES", [ "esint=long", "esint_mpi=MPI_LONG" ])

    ctx.env.append_unique("CXXFLAGS", [ "-std=c++11", "-Wall" ])
    ctx.env.append_unique("CXXFLAGS", ctx.options.cxxflags.split())
    if ctx.options.mode == "release":
        ctx.env.append_unique("CXXFLAGS", [ "-O3", "-g" ])
    if ctx.options.mode == "devel":
        ctx.env.append_unique("CXXFLAGS", [ "-O2", "-g" ])
    if ctx.options.mode == "debug":
        ctx.env.append_unique("CXXFLAGS", [ "-O0", "-g" ])
    if ctx.options.mode == "profile":
        ctx.env.append_unique("CXXFLAGS", [ "-O3" ])
        ctx.options.without_async = True;

    ctx.env.append_unique("INCLUDES", "src")
    ctx.env.append_unique("DEFINES", [ "__ESMODE__="+ctx.options.mode.upper() ])

    """ Recurse to third party libraries wrappers"""
    recurse(ctx)

    if ctx.options.solver.upper() == "PARDISO" and "HAVE_PARDISO" not in ctx.env.DEFINES_PARDISO:
        ctx.options.solver = "mkl"
        Logs.error("Cannot find PARDISO library. Set --pardiso=PATH_TO_PARDISO to use PARDISO solver")
    ctx.env["DEFINES_SOLVER"] = [ "SOLVER_" + ctx.options.solver.upper() ]

    ctx.msg("Setting compiler to", ctx.options.mpicxx)
    ctx.msg("Setting int width to", ctx.options.intwidth)
    ctx.msg("Setting build mode to", ctx.options.mode)
    ctx.msg("Setting solver to", ctx.options.solver)

    print_available(ctx)

def show(ctx):
    ctx.logger = logging.getLogger('show')
    ctx.logger.handlers = Logs.log_handler()

    ctx.msg("CXX", " ".join(ctx.env.CXX))
    ctx.msg("INFO", " ".join(ctx.env.DEFINES_INFO))
    ctx.msg("DEFINES", " ".join(ctx.env.DEFINES))
    ctx.msg("CXXFLAGS", " ".join(ctx.env.CXXFLAGS))
    print_available(ctx)


def env(ctx):
    print(ctx.env)

fetisources = (
   "src/feti/dataholder.cpp",
   "src/feti/generic/Domain.cpp",
   "src/feti/generic/SparseMatrix.cpp",
   "src/feti/generic/utils.cpp",
   "src/feti/generic/timeeval.cpp",
   "src/feti/generic/FETISystemSolver.cpp",
   "src/feti/specific/cluster.cpp",
   "src/feti/specific/itersolver.cpp",
   "src/feti/specific/cpu/clustercpu.cpp",
   "src/feti/specific/cpu/itersolvercpu.cpp",
   "src/feti/specific/cpu/DenseSolverMKL.cpp",
)

def build(ctx):
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip()
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCOMMIT__=\"{0}\"'.format(commit) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCXX__=\"{0}\"'.format(ctx.env.CXX[0]) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESBUILDPATH__=\"{0}\"'.format(ctx.bldnode.abspath()) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCXXFLAGS__=\"{0}\"'.format(" ".join(ctx.env.CXXFLAGS)) ])

    # dirty hack
    # find better solution by waf
    ctx.env["STLIB_MARKER"] = ["-Wl,-Bstatic,--start-group"]
    ctx.env["SHLIB_MARKER"] = ["-Wl,--end-group", "-Wl,-Bdynamic"]

    if ctx.env["DEFINES_SOLVER"][0] == "SOLVER_MKL":
        feti = fetisources + ("src/feti/specific/cpu/SparseSolverMKL.cpp",)
    if ctx.env["DEFINES_SOLVER"][0] == "SOLVER_PARDISO":
        feti = fetisources + ("src/feti/specific/cpu/SparseSolverPARDISO.cpp",)

    features = "cxx cxxshlib"
    ctx.lib = ctx.shlib
    if ctx.options.static:
        features = "cxx"
        ctx.lib = ctx.stlib

    def build(files, target, use=[]):
        ctx(features=features, source=files,target=target, use=use)
        return [ target ] + use

    checker  = build(ctx.path.ant_glob('src/esinfo/**/*.cpp'), "esinfo", [ "INFO" ])
    checker += build(ctx.path.ant_glob('src/config/**/*.cpp'), "config")
    checker += build(ctx.path.ant_glob('src/basis/**/*.cpp'), "basis")

    mesio = list(checker)
    mesio += build(ctx.path.ant_glob('src/mesh/**/*.cpp'), "mesh")
    mesio += build(ctx.path.ant_glob('src/input/**/*.cpp'), "input")
    mesio += build(ctx.path.ant_glob('src/output/**/*.cpp'), "output", [ "ASYNC" ])
    mesio += build(ctx.path.ant_glob('src/wrappers/catalyst/**/*.cpp'), "wcatalyst", [ "CATALYST" ])
    mesio += build(ctx.path.ant_glob('src/wrappers/hdf5/**/*.cpp'), "whdf5", [ "HDF5" ])
    mesio += build(ctx.path.ant_glob('src/wrappers/metis/**/*.cpp'), "wmetis", [ "METIS" ])
    mesio += build(ctx.path.ant_glob('src/wrappers/parmetis/**/*.cpp'), "wparmetis", [ "PARMETIS" ])
    mesio += build(ctx.path.ant_glob('src/wrappers/scotch/**/*.cpp'), "wscotch", [ "SCOTCH" ])
    mesio += build(ctx.path.ant_glob('src/wrappers/ptscotch/**/*.cpp'), "wptscotch", [ "PTSCOTCH" ])
    mesio += build(ctx.path.ant_glob('src/wrappers/kahip/**/*.cpp'), "wkahip", [ "KAHIP" ])

    espreso  = list(mesio)
    espreso += build(ctx.path.ant_glob('src/physics/**/*.cpp'), "physics")
    espreso += build(ctx.path.ant_glob('src/math/**/*.cpp'), "math")
    espreso += build(ctx.path.ant_glob('src/wrappers/mkl/**/*.cpp'), "wmkl", [ "MKL" ])
    espreso += build(ctx.path.ant_glob('src/wrappers/hypre/**/*.cpp'), "whypre", [ "HYPRE" ])
    espreso += build(ctx.path.ant_glob('src/wrappers/mklpdss/**/*.cpp'), "wmklpdss", [ "MKLPDSS", "MKL" ])
    espreso += build(ctx.path.ant_glob('src/wrappers/pardiso/**/*.cpp'), "wpardiso", [ "PARDISO", "MKL" ])
    espreso += build(ctx.path.ant_glob('src/wrappers/superlu/**/*.cpp'), "wsuperlu", [ "SUPERLU", "MKL" ])
    espreso += build(ctx.path.ant_glob('src/wrappers/wsmp/**/*.cpp'), "wwsmp", [ "WSMP" ])
    espreso += build(ctx.path.ant_glob('src/wrappers/bem/**/*.cpp'), "wbem", [ "BEM" ])
    if ctx.env["HAVE_MATH"]:
        espreso += build(feti, "feti", [ "SOLVER", "PARDISO", "MKL" ])

    ctx.program(source="src/app/ecfchecker.cpp", target="ecfchecker", use=checker)
    ctx.program(source="src/app/mesio.cpp", target="mesio", use=mesio, stlib=ctx.options.stlibs, lib=ctx.options.libs)
    if ctx.env["HAVE_MATH"]:
        ctx.program(source="src/app/espreso.cpp",target="espreso", use=espreso, stlib=ctx.options.stlibs, lib=ctx.options.libs)
        ctx.program(source=["src/api/apitester.cpp", "src/api/apidataprovider.cpp"], target="feti4itester", includes="include", use="feti4i")
        ctx.program(source="src/api/example.cpp", target="feti4iexample", includes="include", use="feti4i")
        ctx.lib(source="src/api/wrapper.cpp",target="feti4i", includes="include", use=espreso, stlib=ctx.options.stlibs, lib=ctx.options.libs)

def options(opt):
    opt.compiler = opt.add_option_group("Compiler options")
    opt.decomposers = opt.add_option_group("Third party graph partition tools")
    opt.math = opt.add_option_group("Third party math libraries")
    opt.solvers = opt.add_option_group("Third party solvers")
    opt.other = opt.add_option_group("Other third party libraries")

    opt.compiler.add_option("--mpicxx",
        action="store",
        type="string",
        default="mpiicpc",
        help="MPI compiler used for building of the library [default: %default]")

    opt.compiler.add_option("--cxx",
        action="store",
        choices=["icpc", "g++"],
        metavar="icpc,g++",
        default="icpc",
        help="C++ compiler (set it in the case of non-standard 'mpicxx' settings) [default: %default]")

    opt.compiler.add_option("--cxxflags",
        action="store",
        type="string",
        default="",
        help="C++ compiler flags (space separated list)")

    opt.compiler.add_option("--stlibs",
        action="store",
        type="string",
        default="",
        help="Additional static libraries")

    opt.compiler.add_option("--libs",
        action="store",
        type="string",
        default="",
        help="Additional dynamic libraries")

    opt.compiler.add_option("--intwidth",
        action="store",
        default="32",
        choices=["32", "64"],
        metavar="32,64",
        help="ESPRESO integer datatype width [default: %default]")

    modes=["release", "devel", "debug", "profile"]
    opt.compiler.add_option("-m", "--mode",
        action="store",
        default="release",
        choices=modes,
        help="ESPRESO build mode: " + ", ".join(modes) + " [default: %default]")

    solvers=["mkl", "pardiso"]
    opt.compiler.add_option("--solver",
        action="store",
        default="mkl",
        choices=solvers,
        help="ESPRESO solver " + ", ".join(solvers) + " [default: %default]")

    opt.compiler.add_option("--static",
        action="store_true",
        default=False,
        help="ESPRESO executable file does not contain dynamic libraries (./waf --static).")

    recurse(opt)

