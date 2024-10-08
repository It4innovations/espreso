
import sys, os, logging, subprocess, types

def configure(ctx):
    ctx.env.with_gui = ctx.options.with_gui
    ctx.env.static = ctx.options.static
    ctx.options.with_cuda = ctx.options.with_cuda or ctx.options.solver.upper() == "CUDA"
    ctx.link_cxx = types.MethodType(link_cxx, ctx)

    ctx.msg("Setting int width to", ctx.options.intwidth)
    ctx.msg("Setting build mode to", ctx.options.mode)
    ctx.msg("Setting solver to", ctx.options.solver)

    """ Set compilers """
    ctx.find_program(ctx.options.mpicxx, var="MPICXX")
    if ctx.env.with_gui:
        ctx.env["COMPILER_CXX"] = ctx.env["CXX"]
        ctx.load("compiler_cxx qt5")
    else:
        ctx.load(ctx.options.cxx)
    ctx.env.CXX = ctx.env.LINK_CXX = ctx.env.MPICXX

    """ Set default compilers flags"""

    ctx.env.append_unique("CXXFLAGS", [ "-fopenmp" ])
    ctx.env.append_unique("LINKFLAGS", [ "-fopenmp" ])

    if ctx.options.intwidth == "32":
        ctx.env.append_unique("DEFINES", [ "esint=int", "esint_mpi=MPI_INT" ])
        ctx.env.append_unique("DEFINES_API", [ "FETI4I_INT_WIDTH=32", "MESIO_INT_WIDTH=32" ])
    if ctx.options.intwidth == "64":
        ctx.env.append_unique("DEFINES", [ "esint=long", "esint_mpi=MPI_LONG" ])
        ctx.env.append_unique("DEFINES_API", [ "FETI4I_INT_WIDTH=64", "MESIO_INT_WIDTH=64" ])

    ctx.env.append_unique("CXXFLAGS", [ "-std=c++11", "-Wall" ])  # -fopenmp
    ctx.env.append_unique("CXXFLAGS", ctx.options.cxxflags.split())
    if ctx.options.mode == "release":
        ctx.env.append_unique("CXXFLAGS", [ "-O3", "-g" ])
    if ctx.options.mode == "devel":
        ctx.env.append_unique("CXXFLAGS", [ "-O2", "-g" ])
    if ctx.options.mode == "debug":
        ctx.env.append_unique("CXXFLAGS", [ "-O0", "-g" ])
    if ctx.options.mode == "profile":
        ctx.env.append_unique("CXXFLAGS", [ "-O3" ])

    ctx.env.append_unique("INCLUDES", "src")
    ctx.env.append_unique("DEFINES", [ "__ESMODE__="+ctx.options.mode.upper() ])

    """ Recurse to third party libraries wrappers"""
    recurse(ctx)

    if ctx.options.solver.upper() == "PARDISO" and "HAVE_PARDISO" not in ctx.env.DEFINES_PARDISO:
        ctx.options.solver = "mkl"
        Logs.error("Cannot find PARDISO library. Set --pardiso=PATH_TO_PARDISO to use PARDISO solver")
    ctx.env["DEFINES_SOLVER"] = [ "SOLVER_" + ctx.options.solver.upper() ]

    if ctx.options.with_nvtx:
        ctx.env.append_unique("CXXFLAGS", [ "-DUSE_NVTX" ]) # NVTX profiling tags

    print_available(ctx)

def build(ctx):
    def is_git_directory(path = '.'):
        return subprocess.call(['git', '-C', path, 'status'], stderr=subprocess.STDOUT, stdout = open(os.devnull, 'w')) == 0

    if is_git_directory():
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip().decode()
    else:
        commit = "unknown"

    ctx.env.append_unique("DEFINES_INFO", [ '__ESCOMMIT__=\"{0}\"'.format(commit) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCXX__=\"{0}\"'.format(ctx.env.CXX[0]) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESBUILDPATH__=\"{0}\"'.format(ctx.bldnode.abspath()) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCXXFLAGS__=\"{0}\"'.format(" ".join(ctx.env.CXXFLAGS)) ])

    # dirty hack
    # find better solution by waf
    ctx.env["STLIB_MARKER"] = ["-Wl,-Bstatic,--start-group"]
    ctx.env.prepend_value("SHLIB_MARKER", "-Wl,--end-group")

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

    if ctx.env["DEFINES_SOLVER"][0] == "SOLVER_MKL":
        feti = fetisources + ("src/feti/specific/cpu/SparseSolverMKL.cpp",)
    if ctx.env["DEFINES_SOLVER"][0] == "SOLVER_PARDISO":
        feti = fetisources + ("src/feti/specific/cpu/SparseSolverPARDISO.cpp",)
    if ctx.env["DEFINES_SOLVER"][0] == "SOLVER_CUDA":
        feti = fetisources + ("src/feti/specific/cpu/SparseSolverMKL.cpp", "src/feti/specific/acc/clusterGPU.cpp", "src/feti/specific/acc/itersolverGPU.cpp",)

    features = "cxx cxxshlib"
    ctx.lib = ctx.shlib
    if ctx.env.static or ctx.options.static:
        features = "cxx"
        ctx.lib = ctx.stlib

    prefix = "nb"
    ctx.checker = []
    ctx.mesio = []
    ctx.espreso = []
    def build_checker(files, target, use=[]):
        ctx(features=features, source=files,target=prefix+target, use=use)
        ctx.checker += [ prefix+target ] + use

    def build_mesio(files, target, use=[]):
        ctx(features=features, source=files,target=prefix+target, use=use)
        ctx.mesio += [ prefix+target ] + use

    def build_espreso(files, target, use=[]):
        ctx(features=features, source=files,target=prefix+target, use=use)
        ctx.espreso += [ prefix+target ] + use

    ctx.build_checker = build_checker
    ctx.build_mesio = build_mesio
    ctx.build_espreso = build_espreso

    ctx.build_checker(ctx.path.ant_glob('src/esinfo/**/*.cpp'), "esinfo", [ "INFO" ])
    ctx.build_checker(ctx.path.ant_glob('src/config/**/*.cpp'), "config")
    ctx.build_checker(ctx.path.ant_glob('src/basis/**/*.cpp'), "basis")
    ctx.build_checker(ctx.path.ant_glob('src/autoopt/**/*.cpp'), "autoopt")
    ctx.build_checker(ctx.path.ant_glob('src/wrappers/mpi/**/*.cpp'), "wmpi")

    ctx.build_mesio(ctx.path.ant_glob('src/mesh/**/*.cpp'), "mesh")
    ctx.build_mesio(ctx.path.ant_glob('src/input/**/*.cpp'), "input")
    ctx.build_mesio(ctx.path.ant_glob('src/output/**/*.cpp'), "output", [ "wpthread" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/pthread/**/*.cpp'), "wpthread", [ "PTHREAD" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/catalyst/**/*.cpp'), "wcatalyst", [ "CATALYST" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/hdf5/**/*.cpp'), "whdf5", [ "HDF5" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/gmsh/**/*.cpp'), "wgmsh", [ "GMSH" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/nglib/**/*.cpp'), "wnglib", [ "NGLIB" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/metis/**/*.cpp'), "wmetis", [ "METIS" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/parmetis/**/*.cpp'), "wparmetis", [ "PARMETIS" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/scotch/**/*.cpp'), "wscotch", [ "SCOTCH" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/ptscotch/**/*.cpp'), "wptscotch", [ "PTSCOTCH" ])
    ctx.build_mesio(ctx.path.ant_glob('src/wrappers/kahip/**/*.cpp'), "wkahip", [ "KAHIP" ])

    ctx.lib(source="src/api/wrapper.mesio.cpp", target="mesioapi", includes="include", use=ctx.checker + ctx.mesio + ["API"], stlib=ctx.options.stlibs, lib=ctx.options.libs)
    ctx.program(source=["src/api/api.mesio.cpp"], target="test.mesio", includes="include", use=ctx.checker + ctx.mesio + ["API", "mesioapi"], stlib=ctx.options.stlibs, lib=ctx.options.libs)

    ctx.build_espreso(ctx.path.ant_glob('src/physics/**/*.cpp'), "physics")
    ctx.build_espreso(ctx.path.ant_glob('src/morphing/**/*.cpp'), "devel")
    ctx.build_espreso(ctx.path.ant_glob('src/math/**/*.cpp'), "math")
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/mkl/**/*.cpp'), "wmkl", [ "MKL" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/cuda/**/*.cpp'), "wcuda", [ "CUDA" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/hypre/**/*.cpp'), "whypre", [ "HYPRE" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/mklpdss/**/*.cpp'), "wmklpdss", [ "MKLPDSS", "MKL" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/pardiso/**/*.cpp'), "wpardiso", [ "PARDISO", "MKL" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/superlu/**/*.cpp'), "wsuperlu", [ "SUPERLU", "MKL" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/wsmp/**/*.cpp'), "wwsmp", [ "WSMP" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/csparse/**/*.cpp'), "wcsparse", [ "CSPARSE" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/bem/**/*.cpp'), "wbem", [ "BEM" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/nvtx/**/*.cpp'), "wnvtx", [ "NVTX" ])
    if ctx.env["HAVE_MATH"]:
        if ctx.env.NVCC:
            ctx.build_espreso(ctx.path.ant_glob('src/feti/specific/acc/**/*.cu'), "cudakernels", [ "CUDA" ])
        ctx.build_espreso(feti, "feti", [ "SOLVER", "PARDISO", "MKL" ])

    ctx.program(source="src/app/ecfchecker.cpp", target="ecfchecker", use=ctx.checker)
    ctx.program(source="src/app/mesio.cpp", target="mesio", use=ctx.checker + ctx.mesio, stlib=ctx.options.stlibs, lib=ctx.options.libs)
    if ctx.env["HAVE_MATH"]:
        ctx.program(source="src/app/espreso.cpp",target="espreso", use=ctx.checker + ctx.mesio + ctx.espreso, stlib=ctx.options.stlibs, lib=ctx.options.libs)

        ctx.lib(source="src/api/wrapper.feti4i.cpp", target="feti4i", includes="include", use=ctx.checker + ctx.mesio + ctx.espreso + ["API"], stlib=ctx.options.stlibs, lib=ctx.options.libs)
        ctx.program(source=["src/api/api.feti4i.cpp", "src/api/api.feti4i.dataprovider.cpp"], target="test.feti4i", includes="include", use=ctx.checker + ctx.mesio + ctx.espreso + ["API", "feti4i"], stlib=ctx.options.stlibs, lib=ctx.options.libs)
        ctx.program(source="src/api/example.feti4i.cpp", target="example.feti4i", includes="include", use=ctx.checker + ctx.mesio + ctx.espreso + ["API", "feti4i"], stlib=ctx.options.stlibs, lib=ctx.options.libs)

    if ctx.env.with_gui:
        ctx.objects(source=ctx.path.ant_glob("**/*.ui"), target="ui")
        ctx(
            features = "qt5 cxx cxxprogram",
            source   = ctx.path.ant_glob(["src/gui/**/*.cpp", "src/app/gui.cpp"]),
            moc      = ctx.path.ant_glob("src/gui/**/*.h"),

            use      = [ "ui" ] + ctx.checker + ctx.mesio + ctx.espreso,
            uselib   = "QT5CORE QT5GUI QT5WIDGETS QT5OPENGL",

            target   = "espresogui"
        )

def options(opt):
    opt.load("compiler_cxx qt5")
    opt.compiler = opt.add_option_group("Compiler options")
    opt.decomposers = opt.add_option_group("Third party graph partition tools")
    opt.math = opt.add_option_group("Third party math libraries")
    opt.solvers = opt.add_option_group("Third party solvers")
    opt.other = opt.add_option_group("Other third party libraries")

    opt.compiler.add_option("--mpicxx",
        action="store",
        type="string",
        metavar="MPICXX",
        default=os.path.basename(os.getenv("MPICXX")) if os.getenv("MPICXX") else "mpic++",
        help="MPI compiler used for building of the library")

    opt.compiler.add_option("--cxx",
        action="store",
        type="string",
        metavar="CXX",
        default=os.path.basename(os.getenv("CXX")) if os.getenv("CXX") else "g++",
        help="C++ compiler")
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

    solvers=["mkl", "pardiso", "cuda"]
    opt.compiler.add_option("--solver",
        action="store",
        default="mkl",
        choices=solvers,
        help="ESPRESO solver " + ", ".join(solvers) + " [default: %default]")

    opt.compiler.add_option("--with-cuda",
        action="store_true",
        default=False,
        help="Build ESPRESO with CUDA support.")

    opt.compiler.add_option("--static",
        action="store_true",
        default=False,
        help="ESPRESO executable file does not contain dynamic libraries.")

    opt.compiler.add_option("--with-gui",
        action="store_true",
        default=False,
        help="Build ESPRESO with GUI (Qt5 is needed).")

    recurse(opt)

def print_available(ctx):
    def _print(msg, err, libs, color="RED"):
        libs = [lib for lib in libs if "HAVE_" + lib.upper() in ctx.env["DEFINES_" + lib.upper()]]
        ctx.start_msg(msg)
        if len(libs) == 0:
            ctx.end_msg(err, color=color)
        else:
            ctx.end_msg("[ " + ", ".join(libs) + " ]", color="BLUE")
        return bool(len(libs))

    ctx.env["HAVE_MESH_GENERATOR"] = _print(
        "Available tools for mesh generation",
        "NOT FOUND",
        [ "gmsh", "nglib" ],
        "YELLOW")

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
        [ "mklpdss", "hypre", "pardiso", "superlu", "wsmp", "csparse" ],
        "YELLOW")
    _print(
        "Available miscellaneous libraries",
        "NOT FOUND",
        [ "pthread", "hdf5", "bem", "catalyst", "cuda" ],
        "YELLOW")

""" Recurse to third party libraries wrappers"""
def recurse(ctx):
    """ MPI library """
    ctx.recurse("src/wrappers/mpi")

    """ Graph partition tools """
    ctx.recurse("src/wrappers/metis")
    ctx.recurse("src/wrappers/parmetis")
    ctx.recurse("src/wrappers/scotch")
    ctx.recurse("src/wrappers/ptscotch")
    ctx.recurse("src/wrappers/kahip")

    """ Math libraries"""
    ctx.recurse("src/wrappers/mkl")
    ctx.recurse("src/wrappers/cuda")

    """ Solvers """
    ctx.recurse("src/wrappers/mklpdss")
    ctx.recurse("src/wrappers/hypre")
    ctx.recurse("src/wrappers/pardiso")
    ctx.recurse("src/wrappers/superlu")
    ctx.recurse("src/wrappers/wsmp")
    ctx.recurse("src/wrappers/csparse")

    """ Other """
    ctx.recurse("src/wrappers/pthread")
    ctx.recurse("src/wrappers/hdf5")
    ctx.recurse("src/wrappers/gmsh")
    ctx.recurse("src/wrappers/nglib")
    ctx.recurse("src/wrappers/bem")
    ctx.recurse("src/wrappers/catalyst")
    ctx.recurse("src/wrappers/nvtx")

from waflib import Logs
from waflib.Build import BuildContext
class ShowConfiguration(BuildContext):
    cmd = "show"
    fun = "show"

class ShowEnv(BuildContext):
    cmd = "env"
    fun = "env"

def show(ctx):
    ctx.logger = logging.getLogger('show')
    ctx.logger.handlers = Logs.log_handler()

    ctx.msg("CXX", " ".join(ctx.env.CXX))
    ctx.msg("DEFINES", " ".join(ctx.env.DEFINES))
    ctx.msg("CXXFLAGS", " ".join(ctx.env.CXXFLAGS))
    print_available(ctx)

def env(ctx):
    print(ctx.env)


def link_cxx(self, *k, **kw):
    includes = []
    libpath = []
    if "root" in kw and os.path.isdir(kw["root"]):
        includes = [ os.path.join(kw["root"], dir) for dir in os.listdir(kw["root"]) if dir.startswith("include") ]
        libpath = [ os.path.join(kw["root"], dir) for dir in os.listdir(kw["root"]) if dir.startswith("lib") ]

    general = dict(uselib_store=kw["name"].upper(), mandatory=False)
    if "mandatory" in kw:
        general["mandatory"] = kw["mandatory"]
    if "use" in kw:
        general["use"] = kw["use"]

    header = dict()
    self.env.stash()
    if "header_name" in kw:
        header = dict(header_name=kw["header_name"], define_name="", defines=["HAVE_" + kw["name"].upper()], includes=includes)
        header.update(general)
        header["msg"] = "Checking for '{0}' header".format(kw["name"])
        if not self.check_cxx(**header):
            self.env.revert()
            return False

        if "fragment" in kw:
            test = dict(execute=True)
            test.update(header)
            inc = [ "#include <{0}>\n".format(h) for h in kw["header_name"].split() ]
            test["fragment"] = "{0}int main(int argc, char** argv) {{ {1} }}".format("".join(inc), kw["fragment"])
            test["msg"] = "Checking for '{0}' settings".format(kw["name"])
            if not self.check_cxx(**test):
                self.env.revert()
                return False

    if "libs" in kw:
        libs = dict(stlib=kw["libs"], libpath=libpath, msg="Checking for '{0}' library".format(kw["name"]))
        libs.update(general)
        libs.update(header)
        if not self.options.static or not self.check_cxx(**libs):
            libs["lib"] = libs["stlib"]
            libs.pop("stlib")
            libs["msg"] = "Checking for '{0}' library".format(kw["name"])
            if not self.check_cxx(**libs):
                self.env.revert()
                return False

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

from waflib.TaskGen import feature, before_method, after_method
@feature('qt5')
@after_method('process_source')
@before_method('apply_incpaths')
def add_includes_paths(self):
   incs = set(self.to_list(getattr(self, 'includes', '')))
   for x in self.compiled_tasks:
       incs.add(x.inputs[0].parent.path_from(self.path))
   self.includes = sorted(incs)

from waflib import TaskGen,Task
from waflib.Tools import c_preproc

class cuda(Task.Task):
    run_str = '${NVCC} ${NVCCFLAGS} ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${CPPPATH_ST:INCPATHS} ${DEFINES_ST:DEFINES} ${CXX_SRC_F}${SRC} ${CXX_TGT_F} ${TGT}'
    color   = 'GREEN'
    ext_in  = ['.h']
    vars    = ['CCDEPS']
    scan    = c_preproc.scan
    shell   = False

@TaskGen.extension(".cu", ".cuda")
def cxx_hook(self, node):
    return self.create_compiled_task("cuda", node)

