
VERSION = 15

import commands
import sys
import os

sys.path.append(os.path.abspath("src/python"))
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

solvers = [ "MKL", "PARDISO", "CUDA", "CUDA_7", "MIC", "MUMPS" ]

espreso_attributes = [
    ("CHECK_ENV", "Set to 1, if you want to test the build configuration.", "choice", [ "0", "1" ]),
    ("INT_WIDTH", "ESPRESO integer datatype width.", "choice", [ "32", "64" ]),
    ("LIBTYPE", "ESPRESO is built to libraries of specified type.", "choice", [ "SHARED", "STATIC" ]),
    ("SOLVER", "ESPRESO internal solver. Default: MKL", "choice", solvers),
    ("DEBUG", "Debug information.", "choice", [ "0", "1" ]),
    ("BUILD_TOOLS", "ESPRESO try to compile external tools. If the compilation is not successful set this attribute to 0 and build tools manually.", "choice", [ "0", "1" ]),
    ("METISLIB", "Name of METIS library.", "string", "name"),
]

from waflib.Configure import conf
@conf
def check_header(self, header):
    self.headers += [ header ]

@conf
def check_lib(self, library, dependencies=[]):
    self.libs[library] = [ library ] + dependencies

@conf
def check_stlib(self, library, dependencies=[]):
    self.stlibs[library] = [ library ] + dependencies


def configure(ctx):
    ctx.headers = []
    ctx.libs = {}
    ctx.stlibs = {}

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

    if ctx.env.LIBTYPE == "SHARED":
        ctx.env.append_unique("CXXFLAGS", "-fPIC")

    if ctx.env.BUILD_TOOLS == "1":
        ctx.recurse("tools")

    for header in [ "mpi.h", "mkl.h", "cilk/cilk.h", "omp.h", "tbb/mutex.h" ]:
        ctx.check_header(header)

    # recurse to basic parts
    ctx.recurse("src/config")
    ctx.recurse("src/basis")
    ctx.recurse("src/mesh")
    ctx.recurse("src/input")
    ctx.recurse("src/output")

    # recurse to ESPRESO solver
    ctx.setenv("solver", ctx.env.derive());
    append_solver_attributes(ctx, compiler_attributes)
    ctx.recurse("src/solver")
    ctx.recurse("src/assembler")
    ctx.recurse("src/app")

    check_environment(ctx)

def build(ctx):

<<<<<<< HEAD
    ctx(
        export_includes = "src/include src/config src/basis src/mesh src/input src/output tools/bem4i/src src/assembler src/solver",
        name            = "espreso_includes"
    )
=======
    try:
        ctx.check_cxx(header_name="omp.h", mandatory=True)
    except ctx.errors.ConfigurationError:
        ctx.fatal("Install Open MP or try configuration for your cluster.\n"
            "Run './waf --help' for more options.")

def data_types(ctx):
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
        msg         = "Checking for 64-bit integers")

    espreso = ctx.path.abspath() + "/include/espreso.h"

    ctx.check(
        fragment=
            '''
            #include "''' + espreso + '''"
            #include <stdio.h>

            int main() {
                #if ESPRESO_LOCAL_INDICES_WIDTH == 32
                    printf("32");
                #elif ESPRESO_LOCAL_INDICES_WIDTH == 64
                    printf("64");
                #endif
            }
            ''',
        execute     = True,
        define_ret  = True,
        define_name = "eslocal",
        errmsg      = "Incorrect user-supplied value for ESPRESO_LOCAL_INDICES_WIDTH",
        msg         = "Checking ESPRESO_LOCAL_INDICES_WIDTH"
    )
    eslocal = ctx.get_define("eslocal").replace("\"", "")
    ctx.define("eslocal", int(eslocal))
    ctx.env.ESLOCAL = int(eslocal)

    ctx.check(
        fragment=
            '''
            #include "''' + espreso + '''"
            #include <stdio.h>

            int main() {
                #if ESPRESO_DUAL_INDICES_WIDTH == 32
                    printf("32");
                #elif ESPRESO_DUAL_INDICES_WIDTH == 64
                    printf("64");
                #endif
            }
            ''',
        execute     = True,
        define_ret  = True,
        define_name = "esdual",
        errmsg      = "Incorrect user-supplied value for ESPRESO_DUAL_INDICES_WIDTH",
        msg         = "Checking ESPRESO_DUAL_INDICES_WIDTH"
    )
    esdual = ctx.get_define("esdual").replace("\"", "")
    ctx.define("esdual", int(esdual))
    ctx.env.ESDUAL = int(esdual)

    ctx.check(
        fragment=
            '''
            #include "''' + espreso + '''"
            #include <stdio.h>

            int main() {
                #if ESPRESO_GLOBAL_INDICES_WIDTH == 32
                    printf("32");
                #elif ESPRESO_GLOBAL_INDICES_WIDTH == 64
                    printf("64");
                #endif
            }
            ''',
        execute     = True,
        define_ret  = True,
        define_name = "esglobal",
        errmsg      = "Incorrect user-supplied value for ESPRESO_GLOBAL_INDICES_WIDTH",
        msg         = "Checking ESPRESO_GLOBAL_INDICES_WIDTH"
    )
    esglobal = ctx.get_define("esglobal").replace("\"", "")
    ctx.define("esglobal", int(esglobal))
    ctx.env.ESGLOBAL = int(esglobal)
>>>>>>> RF: rename esint->eslocal and eslong->esglobal

    # Waf INCLUDES policy is strange -> use export includes
    ctx.env.append_unique("CXXFLAGS", [ "-I" + include for include in ctx.env.INCLUDES ])
    ctx.env.INCLUDES = []

    ctx.ROOT = ctx.path.abspath()

    if ctx.env.BUILD_TOOLS == "1":
        ctx.recurse("tools")
    ctx.add_group()

    ctx.recurse("src/basis")
    ctx.recurse("src/config")
    ctx.recurse("src/mesh")
    ctx.recurse("src/input")
    ctx.recurse("src/output")

    ctx.env = ctx.all_envs["solver"]
    ctx.recurse("src/solver")
    ctx.recurse("src/assembler")
    ctx.recurse("src/app")


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

    system = opt.add_option_group("Systems")
    system.add_option(
        "--cray",
        action="store_true",
        default=False,
        help="Compile for Cray"
    )

