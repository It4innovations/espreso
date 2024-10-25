
import sys, os, logging, subprocess, types, copy

def configure(ctx):
    """ Set compilers """
    if ctx.env.with_gui:
        ctx.env["COMPILER_CXX"] = ctx.env["CXX"]
        ctx.load("compiler_cxx qt5")
    else:
        ctx.load("compiler_cxx")

    detect_intel_oneapi_compiler(ctx)
    detect_nvcpp_compiler(ctx)

    if ctx.options.static:
        ctx.env.append_unique("LIB", [ "dl" ])

    """ Set default compilers flags"""
    if ctx.env.COMPILER_CXX == "g++":
        ctx.env.append_unique("CXXFLAGS", [ "-fopenmp", "-Wno-psabi" ])
        ctx.env.append_unique("LINKFLAGS", [ "-fopenmp" ])
    if ctx.env.COMPILER_CXX == "icpc":
        ctx.env.append_unique("CXXFLAGS", [ "-qopenmp", "-diag-disable=10441" ])
        ctx.env.append_unique("LINKFLAGS", [ "-qopenmp", "-diag-disable=10441" ])
    if ctx.env.COMPILER_CXX == "clang++":
        ctx.env.append_unique("CXXFLAGS", [ "-fopenmp" ])
        ctx.env.append_unique("LINKFLAGS", [ "-fopenmp" ])
    if ctx.env.COMPILER_CXX == "icpx":
        ctx.env.append_unique("CXXFLAGS", [ "-qopenmp" ])
        ctx.env.append_unique("LINKFLAGS", [ "-qopenmp" ])
    if ctx.env.COMPILER_CXX == "nvc++":
        ctx.env.append_unique("CXXFLAGS", [ "-fopenmp" ])
        ctx.env.append_unique("LINKFLAGS", [ "-fopenmp" ])

    if ctx.options.flavor == "fujitsu":
        ctx.env.append_unique("CXXFLAGS" , [ "-Kopenmp", "-SSL2" ])
        ctx.env.append_unique("LINKFLAGS", [ "-Kopenmp", "-SSL2" ])

    if ctx.options.intwidth == "32":
        ctx.env.append_unique("DEFINES", [ "esint=int" ])
        ctx.env.append_unique("DEFINES_API", [ "FETI4I_INT_WIDTH=32", "MESIO_INT_WIDTH=32" ])
    if ctx.options.intwidth == "64":
        ctx.env.append_unique("DEFINES", [ "esint=long" ])
        ctx.env.append_unique("DEFINES_API", [ "FETI4I_INT_WIDTH=64", "MESIO_INT_WIDTH=64" ])

    if ctx.options.simd_off:
        ctx.env.append_unique("DEFINES", [ "SIMD_OFF" ])

    ctx.env.append_unique("CXXFLAGS", [ "-std=c++17", "-Wall" ])
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

    ctx.env.intwidth = ctx.options.intwidth
    ctx.env.mode = ctx.options.mode
    settings(ctx)

def build(ctx):
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCOMMIT__=\"{0}\"'.format(get_commit()) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCXX__=\"{0}\"'.format(ctx.env.CXX[0]) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESBUILDPATH__=\"{0}\"'.format(ctx.bldnode.abspath()) ])
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCXXFLAGS__=\"{0}\"'.format(" ".join(ctx.env.CXXFLAGS)) ])

    # dirty hack
    # find better solution by waf
    ctx.env["STLIB_MARKER"] = ["-Wl,-Bstatic,--start-group"]
    ctx.env.prepend_value("SHLIB_MARKER", "-Wl,--end-group")

    features = "cxx cxxshlib"
    ctx.lib = ctx.shlib
    if ctx.env.static or ctx.options.static:
        features = "cxx"
        ctx.lib = ctx.stlib

    prefix = "es"
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
    ctx.build_checker(ctx.path.ant_glob('src/wrappers/exprtk/**/*.cpp'), "wexprtk")
    ctx.build_checker(ctx.path.ant_glob('src/wrappers/backward-cpp/**/*.cpp'), "wbackward", [ "BACKWARD" ])
    ctx.build_checker(ctx.path.ant_glob('src/wrappers/mpi/**/*.cpp'), "wmpi")
    ctx.build_checker(ctx.path.ant_glob('src/wrappers/papi/**/*.cpp'), "wpapi", [ "PAPI" ])

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

    ctx.build_espreso(ctx.path.ant_glob('src/analysis/**/*.cpp'), "analysis")
    ctx.build_espreso(ctx.path.ant_glob('src/morphing/**/*.cpp'), "devel")
    ctx.build_espreso(ctx.path.ant_glob('src/math/**/*.cpp'), "math", [ "BLAS", "LAPACK", "MKL", "SUITESPARSE" ])
    ctx.build_espreso(ctx.path.ant_glob('src/autoopt/**/*.cpp'), "autoopt")
    ctx.build_espreso(ctx.path.ant_glob('src/feti/**/*.cpp'), "feti")
    ctx.build_espreso(ctx.path.ant_glob('src/gpu/**/*.cpp'), "gpu", [ "CUDA", "ROCM", "ONEAPI" ])
#    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/simd/**/*.cpp'), "simd")
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/blas/**/*.cpp'), "wblas", [ "BLAS", "MKL" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/lapack/**/*.cpp'), "wlapack", [ "LAPACK", "MKL" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/mkl/**/*.cpp'), "wmkl", [ "MKL" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/precice/**/*.cpp'), "wprecice", [ "PRECICE" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/precice-dummy/**/*.cpp'), "wprecicedummy", [ "PRECICE" ])
    if ctx.env.HAVE_CUDA:
        ctx.build_espreso(ctx.path.ant_glob('src/wrappers/cuda/**/*.(cu|cpp)'), "wcuda", [ "CUDA" ])
#     ctx.build_espreso(ctx.path.ant_glob('src/wrappers/hypre/**/*.cpp'), "whypre", [ "HYPRE" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/mklpdss/**/*.cpp'), "wmklpdss", [ "MKLPDSS", "MKL" ])
#     ctx.build_espreso(ctx.path.ant_glob('src/wrappers/pardiso/**/*.cpp'), "wpardiso", [ "PARDISO" ])
#     ctx.build_espreso(ctx.path.ant_glob('src/wrappers/superlu/**/*.cpp'), "wsuperlu", [ "SUPERLU" ])
#     ctx.build_espreso(ctx.path.ant_glob('src/wrappers/wsmp/**/*.cpp'), "wwsmp", [ "WSMP" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/suitesparse/**/*.cpp'), "wsuitesparse", [ "SUITESPARSE", "MKL" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/openlb/**/*.cpp'), "wopenlb", [ "OPENLB" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/bem/**/*.cpp'), "wbem",)
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/nvtx/**/*.cpp'), "wnvtx", [ "NVTX" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/oneapi/**/*.cpp'), "woneapi", [ "ONEAPI" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/rocm/**/*.cpp'), "wrocm", [ "ROCM" ])

    ctx.program(source="src/app/ecfchecker.cpp", target="ecfchecker", use=ctx.checker, stlib=ctx.options.stlibs, lib=ctx.options.libs)
    ctx.program(source="src/app/mesio.cpp", target="mesio", use=ctx.checker + ctx.mesio, stlib=ctx.options.stlibs, lib=ctx.options.libs)
    ctx.program(source="src/app/espreso.cpp",target="espreso", rpath=["$ORIGIN/"], use=ctx.checker + ctx.mesio + ctx.espreso, stlib=ctx.options.stlibs, lib=ctx.options.libs)
    ctx.program(source="src/app/dummycoupler.cpp",target="dummycoupler", rpath=["$ORIGIN/"], use=ctx.checker + ctx.mesio + ctx.espreso, stlib=ctx.options.stlibs, lib=ctx.options.libs)

#         ctx.lib(source="src/api/wrapper.feti4i.cpp", target="feti4i", includes="include", use=ctx.checker + ctx.mesio + ctx.espreso + ["API"], stlib=ctx.options.stlibs, lib=ctx.options.libs)
#         ctx.program(source=["src/api/api.feti4i.cpp", "src/api/api.feti4i.dataprovider.cpp"], target="test.feti4i", includes="include", use=ctx.checker + ctx.mesio + ctx.espreso + ["API", "feti4i"], stlib=ctx.options.stlibs, lib=ctx.options.libs)
#         ctx.program(source="src/api/example.feti4i.cpp", target="example.feti4i", includes="include", use=ctx.checker + ctx.mesio + ctx.espreso + ["API", "feti4i"], stlib=ctx.options.stlibs, lib=ctx.options.libs)

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
    opt.load("compiler_cxx")
    opt.compiler = opt.add_option_group("Compiler options")
    opt.decomposers = opt.add_option_group("Third party graph partition tools")
    opt.math = opt.add_option_group("Third party math libraries")
    opt.solvers = opt.add_option_group("Third party sparse solvers")
    opt.other = opt.add_option_group("Other third party libraries")

    opt.compiler.add_option("--cxx",
        action="store",
        type="string",
        metavar="CXX",
        default=os.getenv("CXX") or [ "mpic++" ],
        help="Command to be used for compilation.")

    opt.compiler.add_option("--cxxflags",
        action="store",
        type="string",
        metavar="CXXFLAGS",
        default="",
        help="C++ compiler flags (space separated list)")

    flavor=["gnu", "intel", "fujitsu"]
    opt.compiler.add_option("--flavor",
        action="store",
        default="gnu",
        choices=flavor,
        help="C++ compiler flavor: " + ", ".join(flavor) + " [default: %default]")

    opt.compiler.add_option("--stlibs",
        action="store",
        type="string",
        default="",
        help="Additional static libraries")

    opt.compiler.add_option("--libs",
        action="store",
        type="string",
        default=os.getenv("LIBRARIES"),
        help="Additional dynamic libraries")

    opt.compiler.add_option("--intwidth",
        action="store",
        default=os.getenv("ES_INT_WIDTH") or "32",
        choices=["32", "64"],
        metavar="$ES_INT_WIDTH",
        help="ESPRESO integer datatype width [default: %default]")

    opt.compiler.add_option("--simd-off",
        action="store_true",
        default=False,
        help="Build ESPRESO without SIMD version of assembler.")

    modes=["release", "devel", "debug", "profile"]
    opt.compiler.add_option("-m", "--mode",
        action="store",
        default="release",
        choices=modes,
        help="ESPRESO build mode: " + ", ".join(modes) + " [default: %default]")

    opt.compiler.add_option("--static",
        action="store_true",
        default=False,
        help="ESPRESO executable file does not contain dynamic libraries.")

    opt.compiler.add_option("--with-gui",
        action="store_true",
        default=False,
        help="Build ESPRESO with GUI (Qt5 is needed).")

    opt.other.add_option("--use-cusparse-legacy",
        action="store_true",
        default=os.getenv("ESPRESO_USE_CUSPARSE_LEGACY"),
        help="Use legacy cusparse API. For CUDA < 12 only")

    opt.other.add_option("--rank-to-gpu-map",
        action="store",
        type="string",
        default=os.getenv("ESPRESO_RANK_TO_GPU_MAP") or "0",
        help="Map from node-local MPI rank to GPU index, if the process can see multiple GPUs")

    recurse(opt)

def settings(ctx):
    libraries = { x.lstrip("HAVE_").split("=")[0].lower() for x in ctx.env if x.startswith("HAVE_") }

    def libsmsg(msg, libs):
        libs = [ lib for lib in libs if lib in libraries ]
        ctx.start_msg(msg)
        ctx.end_msg("[ " + ", ".join(libs) + " ]", color="BLUE")

    ctx.msg("                                CXXFLAGS", ctx.env.CXXFLAGS)
    ctx.msg("                               LINKFLAGS", ctx.env.LINKFLAGS)
    if ctx.env.HAVE_CUDA:
        ctx.msg("                                    NVCC", ctx.env.NVCC)
        ctx.msg("                               NVCCFLAGS", ctx.env.NVCCFLAGS)
    ctx.msg("                               int width", ctx.env.intwidth)
    ctx.msg("                                    mode", ctx.env.mode)
    ctx.msg("                                    SIMD", ctx.env.SIMD)
    ctx.msg("                               ASSEMBLER", ctx.env.SIMD_ASSEMBLER)
    libsmsg("                         mesh generators", [ "gmsh", "nglib" ])
    libsmsg("                      graph partitioners", [ "metis", "scotch", "kahip"])
    libsmsg("          distributed graph partitioners", [ "parmetis", "ptscotch" ])
    libsmsg("                          BLAS libraries", [ "mkl", "blas" ])
    libsmsg("                        SpBLAS libraries", [ "mkl", "suitesparse" ])
    libsmsg("                        LAPACK libraries", [ "mkl", "lapack" ])
    libsmsg("                          sparse solvers", [ "mkl", "suitesparse" ])
    libsmsg("              distributed sparse solvers", [ "mklpdss", "hypre", "superlu", "wsmp" ])
    libsmsg("                           GPU libraries", [ "cuda", "rocm", "oneapi" ])
    libsmsg("                                coupling", [ "precice" ])

""" Recurse to third party libraries wrappers"""
def recurse(ctx):
    """ Accelerators """
    ctx.recurse("src/wrappers/rocm")
    ctx.recurse("src/wrappers/cuda")
    ctx.recurse("src/wrappers/oneapi")

    """ Graph partition tools """
    ctx.recurse("src/wrappers/metis")
    ctx.recurse("src/wrappers/parmetis")
    ctx.recurse("src/wrappers/scotch")
    ctx.recurse("src/wrappers/ptscotch")
    ctx.recurse("src/wrappers/kahip")

    """ Math libraries"""
    ctx.recurse("src/wrappers/blas")
    ctx.recurse("src/wrappers/lapack")
    ctx.recurse("src/wrappers/pardiso")
    ctx.recurse("src/wrappers/mkl")
    ctx.recurse("src/wrappers/suitesparse")

    """ Solvers """
    ctx.recurse("src/wrappers/mklpdss")
    ctx.recurse("src/wrappers/hypre")
    ctx.recurse("src/wrappers/superlu")
    ctx.recurse("src/wrappers/wsmp")
    ctx.recurse("src/wrappers/openlb")

    """ Other """
    ctx.recurse("src/wrappers/simd")
    ctx.recurse("src/wrappers/exprtk")
    ctx.recurse("src/wrappers/backward-cpp")
    ctx.recurse("src/wrappers/pthread")
    ctx.recurse("src/wrappers/hdf5")
    # ctx.recurse("src/wrappers/gmsh")
    # ctx.recurse("src/wrappers/nglib")
    # ctx.recurse("src/wrappers/bem")
    # ctx.recurse("src/wrappers/catalyst")
    # ctx.recurse("src/wrappers/nvtx")
    ctx.recurse("src/wrappers/papi")
    ctx.recurse("src/wrappers/precice")

def get_commit():
    def is_git_directory(path = '.'):
        return subprocess.call(['git', '-C', path, 'status'], stderr=subprocess.STDOUT, stdout = open(os.devnull, 'w')) == 0

    if is_git_directory():
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip().decode()
        changed = subprocess.call(['git', 'diff', "--quiet", 'HEAD'], stderr=subprocess.STDOUT, stdout = open(os.devnull, 'w')) == 1
        if changed:
            return commit + "+"
        else:
            return commit
    else:
        return "unknown"

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
    settings(ctx)

def env(ctx):
    print(ctx.env)

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

def detect_intel_oneapi_compiler(ctx):
    ec = os.system("which " + ctx.env.CXX[0] + " > /dev/null 2>/dev/null")
    if ec != 0: return
    ec = os.system(ctx.env.CXX[0] + " --version 2>/dev/null | grep -q \"oneAPI DPC++/C++ Compiler\"")
    if ec != 0: return
    ctx.env.COMPILER_CXX = "icpx"
    
def detect_nvcpp_compiler(ctx):
    ec = os.system(ctx.env.CXX[0] + " --version | grep -q nvc++")
    if ec != 0: return
    ctx.env.COMPILER_CXX = "nvc++"
