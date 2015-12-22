
import os

def check_libraries(ctx):
    def add(map, *libraries):
        for library in libraries:
            map[library] = [ library ]

    libraries = {}
    stlibraries = {}

    add(libraries, "pthread")
    libraries["mkl"] = [ "mkl_intel_lp64", "mkl_core", "mkl_intel_thread" ]

    if ctx.env.SOLVER == "PARDISO":
        libraries["pardiso500-INTEL120-X86-64"] = [ "pardiso500-INTEL120-X86-64", "mkl_intel_lp64", "mkl_core", "mkl_intel_thread", "ifcore" ]
        add(stlibraries, "ifcore")

    if ctx.env.SOLVER == "MIC":
        add(libraries, "imf", "intlc", "svml", "irng")

    if ctx.env.SOLVER == "CUDA":
        add(libraries, "cudart", "cublas")

    if ctx.env.SOLVER == "MUMPS":
        add(libraries, "ifcore")

    for name, library in libraries.items():
        check_library(ctx, "LIB", name, library)

    for name, library in stlibraries.items():
        check_library(ctx, "STLIB", name, library)

def check_headers(ctx):
    headers = [ "mkl.h", "cilk/cilk.h", "omp.h" ]

    for header in headers:
        check_header(ctx, header)

def read_configuration(ctx, espreso_attributes, solvers, compilers, compiler_attributes):
    def read_config(config):
        for line in config:
            parts = line.split('#', 1)[0].partition('=')
            if parts[1] == '=':
                ctx.env[parts[0].strip()] = parts[2].split()

    def read_attribute(attribute, type):
        if getattr(ctx.options, attribute):
            ctx.env[attribute] = getattr(ctx.options, attribute).split()
        if type == "choice":
            ctx.env[attribute] = ctx.env[attribute][0]

    def print_attribute(attribute, type, value):
        if type == "string":
            ctx.msg("Settings " + attribute + " to ", " ".join(value))
        else:
            ctx.msg("Settings " + attribute + " to ", value)

    # Load default configuration
    defaultConfig = open("build.config.default", "r")
    read_config(open("build.config.default", "r"))

    # Load user specific configuration
    if os.path.isfile("build.config"):
        read_config(open("build.config", "r"))

    # Load configuration specified while the project configuration
    for attribute, description, type, value in espreso_attributes + compilers + compiler_attributes:
        read_attribute(attribute, type)

    for attribute, description, type, value in compiler_attributes:
        read_attribute("SOLVER::" + attribute, type)

    # Rewrite solver attributes by attributes from configuration
    for attribute, description, type, value in compiler_attributes:
        if not ctx.env["SOLVER::" + attribute]:
            ctx.env["SOLVER::" + attribute] = ctx.env[ctx.env.SOLVER + "::" + attribute]

    ctx.env.LINK_CXX = ctx.env.CXX

    # Print final configuration
    ctx.find_program(ctx.env.CXX)
    ctx.find_program(ctx.env.CC)
    ctx.find_program(ctx.env.FC)

    for attribute, description, type, value in espreso_attributes:
        print_attribute(attribute, type, ctx.env[attribute])

    for attribute, description, type, value in compiler_attributes:
        print_attribute(attribute, type, ctx.env[attribute])

    for attribute, description, type, value in compiler_attributes:
        print_attribute("SOLVER::" + attribute, type, ctx.env["SOLVER::" + attribute])

def set_datatypes(ctx):
    ctx.env.INT_WIDTH = int(ctx.env.INT_WIDTH)
    if ctx.env.INT_WIDTH == 32:
        ctx.env.append_unique("DEFINES", [ "eslocal=int", "MKL_INT=int", "esglobal=int", "esglobal_mpi=MPI_INT", "BLAS_INT" ])
    elif ctx.env.INT_WIDTH == 64:
        ctx.env.append_unique("DEFINES", [ "eslocal=long", "MKL_INT=long", "esglobal=long", "esglobal_mpi=MPI_LONG", "BLAS_LONG" ])
    else:
        ctx.fatal("ESPRESO supports only INT_WIDTH = {32, 64}.")

def append_solver_attributes(ctx):
    if ctx.env.SOLVER == "MIC" or ctx.env.SOLVER == "CUDA":
        ctx.env.append_unique("DEFINES", ctx.env.SOLVER)

    for attribute in ctx.env.table:
        if attribute.find("SOLVER::") != -1:
            ctx.env.append_unique(attribute.split("::")[1], ctx.env[attribute])

def check_environment(ctx):
    # create new environment and remove all libraries
    ctx.setenv("checker", ctx.env.derive())
    ctx.find_program("cmake")
    ctx.env.LIB = []
    ctx.env.STLIB = []

    ret = ctx.check(
        fragment=
            '''
            #include "mpi.h"
            #include <iostream>
            int main() {
            #ifdef __ICC
                std::cout << __ICC;
                return 0;
            #endif
            #ifdef __INTEL_COMPILER
                std::cout <<  __INTEL_COMPILER;
                return 0;
            #endif
                return 0;
            }
            ''',
        execute     = True,
        define_ret  = True,
        errmsg      = " ".join(ctx.env.CXX) + " is not MPI/Intel compiler",
        msg         = "Checking for MPI/Intel compiler",
        okmsg       = "pass")

    ctx.env.ICPC_VERSION = int(ret[:2])
    ctx.msg("Checking for Intel compiler version", ctx.env.ICPC_VERSION)

    if ctx.env.CHECK_ENV == "0":
        return

    ctx.check_cxx(
        fragment=
            '''
            #include <iostream>
            int main() {
                long x;
                if (sizeof(x) == 8) {
                    std::cout << "sizeof(x)";
                }
                return 0;
            }
            ''',
        execute     = True,
        mandatory   = False,
        errmsg      = "WARNING: ESPRESO with your compiler supports only 32-bit integers",
        okmsg       = "supported",
        msg         = "Checking for 64-bit integers")

    check_headers(ctx)
    check_libraries(ctx)

def check_library(ctx, type, name, library):
    env = ctx.env
    ctx.setenv("tmp", ctx.env.derive());
    ctx.env.append_unique("LIB", library)
    ctx.check(
        msg     = "Checking for library '{0}'".format(name),
        errmsg  = "not found - add path to library '{0}' to {1}PATH".format(name, type),
        okmsg   = "found"
    )
    ctx.env = env

def check_header(ctx, header):
    ctx.check_cxx(
        fragment  = "#include \"{0}\"\nint main() {{ return 0; }}".format(header),
        msg       = "Checking for header '{0}'".format(header),
        errmsg    = "not found - add path to header '{0}' to INCLUDES".format(header),
        okmsg     = "found"
    )
    ctx.env[type] = []

