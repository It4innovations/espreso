
VERSION = 15

import commands
import sys
import os

sys.path.append(os.path.abspath("./python"))
from wafutils import *

# Each attribute has this structure: ( "attribute", "description", "data type", "choices")

compilers = [
    ("CXX", "Intel MPI/C++ compiler used to build ESPRESO.", "string", "compiler"),
    ("CC", "C compiler for build tools used by ESPRESO.", "string", "compiler"),
    ("FC", "Fortran compiler for build tools used by ESPRESO.", "string", "compiler")
]

compiler_attributes = [
    ("CXXFLAGS", "List of compilation flags for ESPRESO files.", "string", "flags"),
    ("LINKFLAGS", "List of flags for linking ESPRESO.", "string", "flags"),
    ("INCLUDES", "Include paths.", "string", "paths"),
    ("LIBPATH", "List of search path for shared libraries.", "string", "paths"),
    ("STLIBPATH", "List of search path for static libraries.", "string", "paths")
]

solvers = [ "MKL", "PARDISO", "CUDA", "MIC", "MUMPS" ]

espreso_attributes = [
    ("CHECK_ENV", "Set to 1, if you want to test the build configuration.", "choice", [ "0", "1" ]),
    ("INT_WIDTH", "ESPRESO integer datatype width.", "choice", [ "32", "64" ]),
    ("LIBTYPE", "ESPRESO is built to libraries of specified type.", "choice", [ "SHARED", "STATIC" ]),
    ("SOLVER", "ESPRESO internal solver. Default: MKL", "choice", solvers),
    ("VERBOSE", "Verbosity level.", "choice", [ "0", "1", "2", "3" ]),
    ("DEBUG", "Debug information.", "choice", [ "0", "1" ]),
    ("BUILD_TOOLS", "ESPRESO try to compile external tools. If the compilation is not successful set this attribute to 0 and build tools manually.", "choice", [ "0", "1" ]),
]


def configure(ctx):
    try:
        ctx.load("icpc")
    except ctx.errors.ConfigurationError:
        ctx.fatal("Install Intel compiler or load the appropriate module.")

    read_configuration(ctx, espreso_attributes, solvers, compilers, compiler_attributes)
    set_compiler_defines(ctx)

    ctx.ROOT = ctx.path.abspath()

    if not os.path.exists(ctx.ROOT + "/libs"):
        os.makedirs(ctx.ROOT + "/libs")

    ctx.env.append_unique("LIBPATH", [ ctx.ROOT + "/libs" ])
    ctx.env.append_unique("STLIBPATH", [ ctx.ROOT + "/libs" ])

    if ctx.env.BUILD_TOOLS == "1":
        ctx.recurse("tools")

    # recurse to basic parts
    ctx.recurse("basis")
    ctx.recurse("config")
    ctx.recurse("bem")
    ctx.recurse("mesh")
    ctx.recurse("input")
    ctx.recurse("output")

    # recurse to ESPRESO solver
    ctx.setenv("solver", ctx.env.derive());
    append_solver_attributes(ctx)
    ctx.recurse("solver")
    ctx.recurse("assembler")
    ctx.recurse("app")

    check_environment(ctx)

def build(ctx):
    ctx(
        export_includes = "basis",
        name            = "incl_basis"
    )
    ctx(
        export_includes = "config",
        name            = "incl_config"
    )
    ctx(
        export_includes = "input",
        name            = "incl_input"
    )
    ctx(
        export_includes = "output",
        name            = "incl_output"
    )
    ctx(
        export_includes = "mesh",
        name            = "incl_mesh"
    )
    ctx(
        export_includes = "solver/src",
        name            = "incl_solver"
    )
    ctx(
        export_includes = "assembler",
        name            = "incl_assembler"
    )
    ctx(
        export_includes = "composer",
        name            = "incl_composer"
    )
    ctx(
        export_includes = "bem/src",
        name            = "incl_bem"
    )
    ctx(
        export_includes = "/usr/local/cuda-7.0/include",
        name            = "incl_cuda"
    )
    ctx(
        export_includes = "include",
        name            = "espreso_includes",
        use             = "incl_basis incl_config incl_input incl_output incl_mesh incl_solver incl_bem incl_catalyst incl_composer incl_assembler incl_cuda"
    )

    ctx.ROOT = ctx.path.abspath()
    ctx.LIBRARIES = ctx.ROOT + "/libs"

    if ctx.env.LIBTYPE == "STATIC":
        ctx.lib = ctx.stlib
    else:
        ctx.lib = ctx.shlib

    ctx.recurse("basis")
    ctx.recurse("config")
    if ctx.env.BUILD_TOOLS == "1":
        ctx.recurse("tools")
    ctx.add_group()
    ctx.recurse("bem")
    ctx.recurse("mesh")
    ctx.recurse("input")
    ctx.recurse("output")

    ctx.env = ctx.all_envs["solver"]
    ctx.recurse("solver")
    ctx.recurse("assembler")
    ctx.recurse("app")


def options(opt):
    opt.parser.formatter.max_help_position = 32

    def add_option(group, attribute, description, type, choices):
        if type == "choice":
            group.add_option("--" + attribute,
            action="store",
            type=type,
            choices=choices,
            metavar="{" + ", ".join(choices) + "}",
            help=description)
        else:
            group.add_option("--" + attribute,
            action="store",
            default="",
            type=type,
            metavar=choices,
            help=description)

    for attribute, description, type, choices in espreso_attributes:
        add_option(
            opt.add_option_group("General ESPRESO parameters"),
            attribute, description, type, choices
        )

    for attribute, description, type, choices in compilers:
        add_option(
            opt.add_option_group("Compilers"),
            attribute, description, type, choices
        )

    for attribute, description, type, choices in compiler_attributes:
        add_option(
            opt.add_option_group("Global compiler attributes"),
            attribute, description, type, choices
        )

    desc = (
        "ESPRESO supports several linear solvers for various platforms. "
        "Each solver is built independently to other ESPRESO parts. "
        "Attributes below specify attributes only for chosen SOLVER.")

    for attribute, description, type, choices in compiler_attributes:
        add_option(
            opt.add_option_group("Solver specific compiler attributes", desc),
            "SOLVER" + "::" + attribute,
            description,
            type, choices
        )



