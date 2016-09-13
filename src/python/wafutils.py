
import os
import subprocess

def check_libraries(ctx):
    def add(map, *libraries):
        for library in libraries:
            map[library] = [ library ]

    def check_library(type, name, library):
        env = ctx.env
        ctx.setenv("tmp", ctx.env.derive());
        ctx.env.append_unique("LIB", library)
        ctx.check(
            msg     = "Checking for library '{0}'".format(name),
            errmsg  = "not found - add path to library '{0}' to {1}PATH".format(name, type),
            okmsg   = "found"
        )
        ctx.env = env

    libraries = {}
    stlibraries = {}

    add(libraries, "z")
    add(libraries, "pthread")
    if ctx.env.BUILD_TOOLS == "0":
        if ctx.env.LIBTYPE == "STATIC":
            add(stlibraries, "metis")
        else:
            add(libraries, "metis")
    libraries["mkl"] = [ "mkl_intel_lp64", "mkl_core", "mkl_intel_thread" ]

    if ctx.env.SOLVER == "PARDISO":
        libraries["pardiso500-INTEL120-X86-64"] = [ "pardiso500-INTEL120-X86-64", "mkl_intel_lp64", "mkl_core", "mkl_intel_thread", "ifcore" ]
        add(stlibraries, "ifcore")

    if ctx.env.SOLVER == "MIC":
        add(libraries, "imf", "intlc", "svml", "irng")

    if ctx.env.SOLVER == "CUDA":
        add(libraries, "cudart", "cublas", "cusolver", "cusparse")

    if ctx.env.SOLVER == "MUMPS":
        add(libraries, "ifcore")

    for name, library in libraries.items():
        check_library("LIB", name, library)

    for name, library in stlibraries.items():
        check_library("STLIB", name, library)

def check_headers(ctx):
    for header in ctx.headers:
        ctx.check_cc(
            header_name = header,
            execute     = False,
            msg         = "Checking for header '{0}'".format(header),
            errmsg      = "not found - add path to header '{0}' to INCLUDES".format(header),
            okmsg       = "found"
        )

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
    if ctx.options.cray:
        read_config(open("build.config.cray", "r"))
    else:
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

    for out_attribute in [ "HYPRE::INCLUDE", "HYPRE::LIBPATH", "VTK::INCLUDE", "VTK::LIBPATH", "PARAVIEW::INCLUDE", "PARAVIEW::LIBPATH"]:
        print_attribute(out_attribute, "string", ctx.env[out_attribute])

    for attribute, description, type, value in compiler_attributes:
        print_attribute("SOLVER::" + attribute, type, ctx.env["SOLVER::" + attribute])

def set_compiler_defines(ctx):
    ctx.env.INT_WIDTH = int(ctx.env.INT_WIDTH)
    if ctx.env.INT_WIDTH == 32:
        ctx.env.append_unique("DEFINES", [ "eslocal=int", "MKL_INT=int", "esglobal=int", "esglobal_mpi=MPI_INT", "BLAS_INT", "INT_WIDTH=32" ])
    elif ctx.env.INT_WIDTH == 64:
        ctx.env.append_unique("DEFINES", [ "eslocal=long", "MKL_INT=long", "esglobal=long", "esglobal_mpi=MPI_LONG", "BLAS_LONG", "INT_WIDTH=64" ])
    else:
        ctx.fatal("ESPRESO supports only INT_WIDTH = {32, 64}.")

    ctx.env.append_unique("DEFINES", [ "FETI4I_INT_WIDTH={0}".format(ctx.env.INT_WIDTH) ])

    if ctx.env.DEBUG == "1":
        ctx.env.append_unique("DEFINES", [ "DEBUG" ])

def append_solver_attributes(ctx, attributes):
    if ctx.env.SOLVER == "MIC" or ctx.env.SOLVER == "CUDA":
        ctx.env.append_unique("DEFINES", ctx.env.SOLVER)
    ctx.env.append_unique("DEFINES", [ "SOLVER_{0}".format(ctx.env.SOLVER) ])

    for attribute in attributes:
        ctx.env.append_unique(attribute[0], ctx.env["SOLVER::" + attribute[0]])

def check_environment(ctx):
    if ctx.env.CHECK_ENV == "0":
        return

    # create new environment and remove all libraries
    ctx.setenv("checker", ctx.env.derive())
    ctx.find_program("cmake")
    ctx.env.LIB = []
    ctx.env.STLIB = []

    ret = subprocess.Popen(ctx.env.CXX + ctx.env.CXXFLAGS + [ "--version" ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    version = ret.communicate()[0].split()[2]
    ctx.env.ICPC_VERSION = int(version[:2])
    ctx.msg("Checking for Intel compiler version", ctx.env.ICPC_VERSION)

    check_headers(ctx)
    check_libraries(ctx)


