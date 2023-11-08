
import sys, os, logging, subprocess, types

libs_blas=[ "cblas", "openblas", "sslblas", "mkl" ]
libs_spblas=[ "mkl", "suitesparse" ]
libs_lapack=[ "mkl", "lapacke", "ssllapack" ]
libs_solvers=[ "mkl", "pardiso", "suitesparse" ]

def configure(ctx):
    ctx.env.with_gui = ctx.options.with_gui
    ctx.env.static = ctx.options.static
    ctx.test_dict = types.MethodType(test_dict, ctx)
    ctx.link_cxx = types.MethodType(link_cxx, ctx)
    ctx.env.intwidth = ctx.options.intwidth
    ctx.env.mode = ctx.options.mode

    if ctx.options.cxx == "icpc" or ctx.options.cxx == "icx":
        ctx.options.flavor = "intel"

    """ Set compilers """
    ctx.find_program(ctx.options.mpicxx, var="MPICXX")
    if ctx.env.with_gui:
        ctx.env["COMPILER_CXX"] = ctx.env["CXX"]
        ctx.load("compiler_cxx qt5")
    else:
        ctx.load(ctx.options.cxx)
    ctx.env.CXX = ctx.env.LINK_CXX = ctx.env.MPICXX

    """ Set default compilers flags"""
    if ctx.options.flavor == "gnu":
        ctx.env.append_unique("CXXFLAGS", [ "-fopenmp" ])
        ctx.env.append_unique("LINKFLAGS", [ "-fopenmp" ])
    if ctx.options.flavor == "intel":
        ctx.env.append_unique("CXXFLAGS", [ "-qopenmp", "-diag-disable=10441" ])
        ctx.env.append_unique("LINKFLAGS", [ "-qopenmp", "-diag-disable=10441" ])
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

    ctx.env.append_unique("CXXFLAGS", [ "-std=c++17" ])
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
    if ctx.options.use_spblas == "suitesparse" or ctx.options.use_solver == "suitesparse":
        ctx.options.use_spblas = "suitesparse"
        ctx.options.use_solver = "suitesparse"
    ctx.options.use_blas = set_libs(ctx, libs_blas, ctx.options.use_blas)
    ctx.options.use_spblas = set_libs(ctx, libs_spblas, ctx.options.use_spblas)
    ctx.options.use_lapack = set_libs(ctx, libs_lapack, ctx.options.use_lapack)
    ctx.options.use_solver = set_libs(ctx, libs_solvers, ctx.options.use_solver)
    ctx.env.append_unique("DEFINES_"+ctx.options.use_blas.upper(), "USE_BLAS_"+ctx.options.use_blas.upper())
    ctx.env.append_unique("DEFINES_"+ctx.options.use_spblas.upper(), "USE_SPBLAS_"+ctx.options.use_spblas.upper())
    ctx.env.append_unique("DEFINES_"+ctx.options.use_lapack.upper(), "USE_LAPACK_"+ctx.options.use_lapack.upper())
    ctx.env.append_unique("DEFINES_"+ctx.options.use_solver.upper(), "USE_SOLVER_"+ctx.options.use_solver.upper())
    ctx.env.use_blas = ctx.options.use_blas;
    ctx.env.use_spblas = ctx.options.use_spblas;
    ctx.env.use_lapack = ctx.options.use_lapack;
    ctx.env.use_solver = ctx.options.use_solver;

    if ctx.options.with_nvtx:
        ctx.env.append_unique("CXXFLAGS", [ "-DUSE_NVTX" ]) # NVTX profiling tags

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
    ctx.build_checker(ctx.path.ant_glob('src/wrappers/backward-cpp/**/*.cpp'), "wbackward", [ "BACKWARD-CPP" ])
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
#     ctx.build_espreso(ctx.path.ant_glob('src/physics/**/*.cpp'), "physics")
    ctx.build_espreso(ctx.path.ant_glob('src/morphing/**/*.cpp'), "devel")
    ctx.build_espreso(ctx.path.ant_glob('src/math/**/*.cpp'), "math", [ "MKL", "PARDISO", "SUITESPARSE", "CBLAS", "OPENBLAS", "SSLBLAS" ])
    ctx.build_espreso(ctx.path.ant_glob('src/autoopt/**/*.cpp'), "autoopt")
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/simd/**/*.cpp'), "simd")
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/cblas/**/*.cpp'), "wcblas", [ "CBLAS" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/openblas/**/*.cpp'), "wopenblas", [ "OPENBLAS" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/sslblas/**/*.cpp'), "wsslblas", [ "SSLBLAS" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/lapacke/**/*.cpp'), "wlapacke", [ "LAPACKE" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/ssllapack/**/*.cpp'), "wssllapack", [ "SSLLAPACK" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/mkl/**/*.cpp'), "wmkl", [ "MKL", "PARDISO" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/cuda/**/*.cpp'), "wcuda", [ "CUDA" ])
#     ctx.build_espreso(ctx.path.ant_glob('src/wrappers/hypre/**/*.cpp'), "whypre", [ "HYPRE" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/mklpdss/**/*.cpp'), "wmklpdss", [ "MKLPDSS", "MKL" ])
#     ctx.build_espreso(ctx.path.ant_glob('src/wrappers/pardiso/**/*.cpp'), "wpardiso", [ "PARDISO", "MKL" ])
#     ctx.build_espreso(ctx.path.ant_glob('src/wrappers/superlu/**/*.cpp'), "wsuperlu", [ "SUPERLU", "MKL" ])
#     ctx.build_espreso(ctx.path.ant_glob('src/wrappers/wsmp/**/*.cpp'), "wwsmp", [ "WSMP" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/csparse/**/*.cpp'), "wcsparse", [ "CSPARSE" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/suitesparse/**/*.cpp'), "wsuitesparse", [ "SUITESPARSE" ])
#     ctx.build_espreso(ctx.path.ant_glob('src/wrappers/bem/**/*.cpp'), "wbem", [ "BEM" ])
    ctx.build_espreso(ctx.path.ant_glob('src/wrappers/nvtx/**/*.cpp'), "wnvtx", [ "NVTX" ])
#         if ctx.env.NVCC:
#             ctx.build_espreso(ctx.path.ant_glob('src/feti/specific/acc/**/*.cu'), "cudakernels", [ "CUDA" ])
#         ctx.build_espreso(feti, "feti", [ "SOLVER", "PARDISO", "MKL" ])
    ctx.build_espreso(ctx.path.ant_glob('src/feti/**/*.cpp'), "feti")

    ctx.program(source="src/app/ecfchecker.cpp", target="ecfchecker", use=ctx.checker)
    ctx.program(source="src/app/mesio.cpp", target="mesio", use=ctx.checker + ctx.mesio, stlib=ctx.options.stlibs, lib=ctx.options.libs)
    ctx.program(source="src/app/espreso.cpp",target="espreso", rpath=["$ORIGIN/"], use=ctx.checker + ctx.mesio + ctx.espreso, stlib=ctx.options.stlibs, lib=ctx.options.libs)

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
    opt.load("compiler_cxx qt5")
    opt.compiler = opt.add_option_group("Compiler options")
    opt.decomposers = opt.add_option_group("Third party graph partition tools")
    opt.math = opt.add_option_group("Third party math libraries")
    opt.solvers = opt.add_option_group("Third party sparse solvers")
    opt.other = opt.add_option_group("Other third party libraries")

    opt.math.add_option("--use-blas", action="store", default="mkl", choices=libs_blas, help="Set BLAS library to use: {} [default: %default].".format(libs_blas))
    opt.math.add_option("--use-spblas", action="store", default="mkl", choices=libs_spblas, help="Set SpBLAS library to use: {} [default: %default].".format(libs_spblas))
    opt.math.add_option("--use-lapack", action="store", default="mkl", choices=libs_lapack, help="Set LAPACk library to use: {} [default: %default].".format(libs_lapack))
    opt.solvers.add_option("--use-solver", action="store", default="mkl", choices=libs_solvers, help="Set sparse solver library to use: {} [default: %default].".format(libs_solvers))

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
        default="",
        help="Additional dynamic libraries")

    opt.compiler.add_option("--intwidth",
        action="store",
        default="32",
        choices=["32", "64"],
        metavar="32,64",
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

    recurse(opt)

def set_libs(ctx, libs, default):
    use = default
    for lib in libs:
        if "HAVE_"+lib.upper() in ctx.env["DEFINES_"+lib.upper()]:
            use = lib
            if lib == default:
                return lib
    return use

def settings(ctx):
    print_available(ctx)

    ctx.msg("                                CXXFLAGS", ctx.env.CXXFLAGS)
    ctx.msg("                               LINKFLAGS", ctx.env.LINKFLAGS)
    ctx.msg("                               int width", ctx.env.intwidth)
    ctx.msg("                                    SIMD", ctx.env.SIMD)
    ctx.msg("                                    mode", ctx.env.mode)
    ctx.msg("                                    BLAS", ctx.env.use_blas)
    ctx.msg("                                  SpBLAS", ctx.env.use_spblas)
    ctx.msg("                                  LAPACK", ctx.env.use_lapack)
    ctx.msg("                           sparse solver", ctx.env.use_solver)
    ctx.env.mesio = ctx.env.decomposers or ctx.env.dist_decomposers
    ctx.env.espreso = ctx.env.blas and ctx.env.spblas and ctx.env.solvers
    if ctx.env.mesio and ctx.env.espreso:
        ctx.msg("            MESIO, ESPRESO configuration", "passed");
    else:
        if ctx.env.decomposers or ctx.env.dist_decomposers:
            ctx.msg("                     MESIO configuration", "passed");
        else:
            ctx.msg("                     MESIO configuration", "partialy passed");
        if ctx.env.blas and ctx.env.spblas and ctx.env.solvers:
            ctx.msg("                   ESPRESO configuration", "passed");
        else:
            ctx.msg("                   ESPRESO configuration", "failed");

def print_available(ctx):
    def _print(msg, libs):
        libs = [lib for lib in libs if "HAVE_" + lib.upper() in ctx.env["DEFINES_" + lib.upper()]]
        ctx.start_msg(msg)
        ctx.end_msg("[ " + ", ".join(libs) + " ]", color="BLUE")
        return len(libs)

    ctx.env.generators       = _print("                         mesh generators", [ "gmsh", "nglib" ])
    ctx.env.decomposers      = _print("                      graph partitioners", [ "metis", "scotch", "kahip"])
    ctx.env.dist_decomposers = _print("          distributed graph partitioners", [ "parmetis", "ptscotch" ])
    ctx.env.blas             = _print("                          BLAS libraries", libs_blas)
    ctx.env.spblas           = _print("                        SpBLAS libraries", libs_spblas)
    ctx.env.lapack           = _print("                        LAPACK libraries", libs_lapack)
    ctx.env.solvers          = _print("                          sparse solvers", libs_solvers)
    ctx.env.dist_solvers     = _print("              distributed sparse solvers", [ "mklpdss", "hypre", "superlu", "wsmp" ])
    ctx.env.other            = _print("                 miscellaneous libraries", [ "hdf5", "pthread", "PAPI" ])

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
    ctx.recurse("src/wrappers/cblas")
    ctx.recurse("src/wrappers/openblas")
    ctx.recurse("src/wrappers/sslblas")
    ctx.recurse("src/wrappers/lapacke")
    ctx.recurse("src/wrappers/ssllapack")
    ctx.recurse("src/wrappers/pardiso")
    ctx.recurse("src/wrappers/mkl")
    ctx.recurse("src/wrappers/suitesparse")
    ctx.recurse("src/wrappers/cuda")
    ctx.recurse("src/wrappers/csparse")

    """ Solvers """
    ctx.recurse("src/wrappers/mklpdss")
    ctx.recurse("src/wrappers/hypre")
    ctx.recurse("src/wrappers/superlu")
    ctx.recurse("src/wrappers/wsmp")

    """ Other """
    ctx.recurse("src/wrappers/simd")
    ctx.recurse("src/wrappers/exprtk")
    ctx.recurse("src/wrappers/backward-cpp")
    ctx.recurse("src/wrappers/pthread")
    ctx.recurse("src/wrappers/hdf5")
    ctx.recurse("src/wrappers/gmsh")
    ctx.recurse("src/wrappers/nglib")
    ctx.recurse("src/wrappers/bem")
    ctx.recurse("src/wrappers/catalyst")
    ctx.recurse("src/wrappers/nvtx")
    ctx.recurse("src/wrappers/papi")

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

class GetInfo(BuildContext):
    cmd = "info"
    fun = "info"

def show(ctx):
    ctx.logger = logging.getLogger('show')
    ctx.logger.handlers = Logs.log_handler()
    settings(ctx)

def env(ctx):
    print(ctx.env)

def info(ctx):
    parameters = [
        ("commit", get_commit()),
        ("int_width", ctx.env.intwidth),
        ("mode", ctx.env.mode),
        ("BLAS", ctx.env.use_blas),
        ("SpBLAS", ctx.env.use_spblas),
        ("LAPACK", ctx.env.use_lapack),
        ("solver", ctx.env.use_solver)
        ]
    print("_".join(map(lambda setting: "{}-{}".format(setting[0], setting[1]), parameters)))

def test_dict(self, name):
    test = dict()
    test["msg"] = "Checking for " + name
    test["define_name"] = ""
    test["defines"] = "HAVE_" + name.upper()
    test["name"] = test["uselib_store"] = name.upper()
    return test

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
        header = dict(header_name=kw["header_name"], define_name="", includes=includes)
        header.update(general)
        header["msg"] = "Checking for '{0}' header".format(kw["name"])
        if not self.check_cxx(**header):
            self.env.revert()
            return False

        if "fragment" in kw:
            test = dict(execute=True)
            test.update(header)
            inc = [ "#include <{0}>\n".format(h) for h in kw["header_name"].split() ]
            test["fragment"] = "{0}int main(int argc, char** argv) {{ {1} }}\n".format("".join(inc), kw["fragment"])
            test["msg"] = "Checking for '{0}' settings".format(kw["name"])
            if not self.check_cxx(**test):
                self.env.revert()
                return False

    if "libs" in kw:
        libs = dict(stlib=kw["libs"], libpath=libpath, msg="Checking for '{0}' library".format(kw["name"]))
        libs.update(general)
        libs.update(header)
        libs["msg"] = "Checking for '{0}' library".format(kw["name"])
        if not self.options.static or not self.check_cxx(**libs):
            libs["lib"] = libs["stlib"]
            libs.pop("stlib")
            if not self.check_cxx(**libs):
                self.env.revert()
                return False

    self.env.append_unique("DEFINES_"+kw["name"].upper(), "HAVE_" + kw["name"].upper())
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

