
def options(opt):
    opt.load("icpc")

    opt.add_option("--debug",
       action="store_true",
       default=False,
       help="Compile sequential code.")

    opt.add_option("--mpich",
       action="store_true",
       default=False,
       help="Compile MPI version with mpich.")

    opt.add_option("--mesh",
       action="store_true",
       default=False,
       help="Create application only from mesh.")


def check_environment(ctx):
    try:
        ctx.load("icpc")
    except ctx.errors.ConfigurationError:
        ctx.fatal("Install Intel compiler or try configuration for your cluster.\n"
            "Run './waf --help' for more options.")

    try:
        ctx.find_program("mpic++", var="MPICXX")
    except ctx.errors.ConfigurationError:
        ctx.fatal("mpic++ not found. Install MPI or try configuration for your cluster.\n"
            "Run './waf --help' for more options.")

    try:
        ctx.check_cxx(header_name="mkl.h", mandatory=True)
        ctx.check_cxx(header_name="cilk/cilk.h", mandatory=True)
        ctx.check_cxx(header_name="tbb/tbb.h", mandatory=True)
    except ctx.errors.ConfigurationError:
        ctx.fatal("Install MKL or try configuration for your cluster.\n"
            "Run './waf --help' for more options.")

    try:
        ctx.check_cxx(header_name="omp.h", mandatory=True)
    except ctx.errors.ConfigurationError:
        ctx.fatal("Install Open MP or try configuration for your cluster.\n"
            "Run './waf --help' for more options.")


def configure(ctx):
    check_environment(ctx)

    ctx.setenv("base", ctx.env)
    ctx.env.append_unique("CXXFLAGS", [ "-Wall" ])
    ctx.env.append_unique("LIB", [ "tbb" ])

    if ctx.options.debug:
        ctx.env.append_unique("CXXFLAGS", [ "-g", "-mkl=sequential", "-openmp", "-DSEQUENTIAL" ])
        ctx.env.append_unique("LINKFLAGS", [ "-mkl=sequential", "-openmp" ])
        ctx.env.append_unique("LIB", [ "gomp" ])
    else:
        ctx.env.append_unique("CXXFLAGS", [ "-O2", "-mkl=parallel" ])
        ctx.env.append_unique("LINKFLAGS", [ "-mkl=parallel", "-Wl,-rpath,../libs" ])

    ctx.setenv("mpi", ctx.env)
    if ctx.options.mpich:
        ctx.env.append_unique("CXXFLAGS", [ "-cxx=icpc" ])
        ctx.env.append_unique("LINKFLAGS", [ "-cxx=icpc" ])

    ctx.env.CXX = list(ctx.env.MPICXX)
    ctx.env.LINK_CXX = list(ctx.env.MPICXX)

    ctx.recurse("bem")
    ctx.recurse("mesh")
    ctx.recurse("permoncube")
    ctx.recurse("solver")
    ctx.recurse("app")

def build(ctx):
    ctx(
        export_includes = "../metis/metis-5.1.0/include",
        name            = "metis_includes"
    )

    ctx.ROOT = ctx.path.abspath()

    ctx.recurse("bem")
    ctx.recurse("mesh")

    if ctx.options.mesh:
        return

    ctx.recurse("permoncube")
    ctx.recurse("solver")
    ctx.recurse("app")
