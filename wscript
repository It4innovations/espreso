
def options(opt):
    opt.load("icpc")

    opt.add_option("--debug",
       action="store_true",
       default=False,
       help="Compile sequential code")

    opt.add_option("--mpich",
       action="store_true",
       default=False,
       help="Compile MPI version with mpich")

    opt.add_option("--openmpi",
       action="store_true",
       default=False,
       help="Compile MPI version with openMPI")

def configure(ctx):
    ctx.load("icpc")

    ctx.setenv("base", ctx.env)
    ctx.env.append_unique("CXXFLAGS", [ "-O2"])
    ctx.env.append_unique("LIB", [ "tbb" ])

    if ctx.options.debug:
        ctx.env.append_unique("CXXFLAGS", [ "-g", "-mkl=sequential" ])
        ctx.env.append_unique("LINKFLAGS", [ "-mkl=sequential" ])
    else:
        ctx.env.append_unique("CXXFLAGS", [ "-mkl=parallel" ])
        ctx.env.append_unique("LINKFLAGS", [ "-mkl=parallel" ])

    ctx.setenv("mpi", ctx.env)
    try:
        if ctx.options.mpich or ctx.options.openmpi:
            ctx.find_program("mpic++", var="MPICXX")
            if ctx.options.mpich:
                ctx.env.append_unique("CXXFLAGS", [ "-cxx=icpc" ])
                ctx.env.append_unique("LINKFLAGS", [ "-cxx=icpc" ])
        else:
            ctx.find_program("mpiicc", var="MPICXX")
    except ctx.errors.ConfigurationError:
        ctx.fatal("MPI not found. Try configuration with --mpich or --openmpi")

    ctx.env.CXX = list(ctx.env.MPICXX)
    ctx.env.LINK_CXX = list(ctx.env.MPICXX)

    ctx.recurse("mesh")
    ctx.recurse("permoncube")
    ctx.recurse("solver")
    ctx.recurse("app")


def build(ctx):
    ctx.recurse("mesh")
    ctx.recurse("permoncube")
    ctx.recurse("solver")
    ctx.recurse("app")





