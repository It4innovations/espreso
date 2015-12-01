
import os

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

    ret = ctx.check(
        fragment=
            '''
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
                std::cout << \"0\";
                return 0;
            }
            ''',
        execute     = True,
        define_ret  = True,
        errmsg      = "Incorrect user-supplied value for ESPRESO_LOCAL_INDICES_WIDTH",
        msg         = "Checking icpc compiler support.",
        okmsg       = "supported"
    )
    ctx.env.ICPC_VERSION = int(ret[:2])
    ctx.msg("Checking icpc version", ctx.env.ICPC_VERSION)

def set_default(ctx):
    try:
        if ctx.options.mpich:
            ctx.find_program("mpic++.mpich", var="MPICXX", mandatory=False)
            if not ctx.env.MPICXX:
                ctx.find_program("mpic++", var="MPICXX")
            ctx.find_program("mpicc.mpich", var="MPICC", mandatory=False)
            if not ctx.env.MPICC:
                ctx.find_program("mpicc", var="MPICC")
            ctx.find_program("mpif90.mpich", var="MPIFORT", mandatory=False)
            if not ctx.env.MPIFORT:
                ctx.find_program("mpif90", var="MPIFORT")
        else:
            ctx.find_program("mpic++", var="MPICXX")
            ctx.find_program("mpiicc", var="MPICC")
            ctx.find_program("mpiifort", var="MPIFORT")

    except ctx.errors.ConfigurationError:
        ctx.fatal("mpic++ not found.\nRun './waf --help' for more options.")

    if ctx.options.mpich:
        ctx.env.append_unique("CXXFLAGS", [ "-cxx=icpc" ])
        ctx.env.append_unique("LINKFLAGS", [ "-cxx=icpc" ])
        ctx.options.gfortran = True

    ctx.env.CXX = ctx.env.LINK_CXX = ctx.env.MPICXX;

    ctx.setenv("tmp", ctx.all_envs[""].derive())

    try:
        ctx.check_cxx(header_name="mkl.h", mandatory=True)
        ctx.check_cxx(header_name="cilk/cilk.h", mandatory=True)
        ctx.check_cxx(header_name="tbb/tbb.h", mandatory=True)
    except ctx.errors.ConfigurationError:
        ctx.fatal("Install MKL.\nRun './waf --help' for more options.")

    try:
        ctx.check_cxx(header_name="omp.h", mandatory=True)
    except ctx.errors.ConfigurationError:
        ctx.fatal("Install Open MP.\nRun './waf --help' for more options.")

    ctx.env = ctx.all_envs[""]



def configure_indices_width(ctx):

    espresoHeader = open("include/espreso.h", "r")
    indices = {}

    for line in espresoHeader:
        if line.startswith("#define"):
            values = line.split()
            if len(values) > 2:
                indices[values[1]] = values[2]

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
        msg         = "Checking ESPRESO_LOCAL_INDICES_WIDTH",
        okmsg       = indices["ESPRESO_LOCAL_INDICES_WIDTH"]
    )
    ctx.define("eslocal", int(indices["ESPRESO_LOCAL_INDICES_WIDTH"]))
    ctx.env.ESLOCAL = int(indices["ESPRESO_LOCAL_INDICES_WIDTH"])

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
        msg         = "Checking ESPRESO_DUAL_INDICES_WIDTH",
        okmsg       = indices["ESPRESO_DUAL_INDICES_WIDTH"]
    )
    ctx.define("esdual", int(indices["ESPRESO_DUAL_INDICES_WIDTH"]))
    ctx.env.ESDUAL = int(indices["ESPRESO_DUAL_INDICES_WIDTH"])

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
        msg         = "Checking ESPRESO_GLOBAL_INDICES_WIDTH",
        okmsg       = indices["ESPRESO_GLOBAL_INDICES_WIDTH"]
    )
    ctx.define("esglobal", int(indices["ESPRESO_GLOBAL_INDICES_WIDTH"]))
    ctx.env.ESGLOBAL = int(indices["ESPRESO_GLOBAL_INDICES_WIDTH"])


def check_configuration(version):
    return True
    if not os.path.isfile("build/config.h") or not os.path.isfile("include/espreso.h"):
        return False

    configFile = open("build/config.h", "r")
    espresoFile = open("include/espreso.h", "r")

    config = {}
    espreso = {}

    for line in espresoFile:
        if line.startswith("#define"):
            values = line.split()
            if len(values) > 2:
                espreso[values[1]] = values[2]

    for line in configFile:
        if line.startswith("#define"):
            values = line.split()
            if len(values) > 2:
                config[values[1]] = values[2]

    if config["eslocal"] != espreso["ESPRESO_LOCAL_INDICES_WIDTH"]:
        return False
    if config["esdual"] != espreso["ESPRESO_DUAL_INDICES_WIDTH"]:
        return False
    if config["esglobal"] != espreso["ESPRESO_GLOBAL_INDICES_WIDTH"]:
        return False
    if version != int(config["version"]):
        return False

    # all settings are the same
    return True


