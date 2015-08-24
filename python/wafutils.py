
import os

def set_indices_width(ctx):
    if ctx.env.ESLOCAL == 32:
        ctx.env.append_unique("CXXFLAGS", [ "-Deslocal=int", "-DMKL_INT=int" ])
    if ctx.env.ESLOCAL == 64:
        ctx.env.append_unique("CXXFLAGS", [ "-Deslocal=long", "-DMKL_INT=long" ])

    if ctx.env.ESDUAL == 32:
        ctx.env.append_unique("CXXFLAGS", [ "-Desdual=int" ])
    if ctx.env.ESDUAL == 64:
        ctx.env.append_unique("CXXFLAGS", [ "-Desdual=long" ])

    if ctx.env.ESGLOBAL == 32:
        ctx.env.append_unique("CXXFLAGS", [ "-Desglobal=int", "-Desglobal_mpi=MPI_INT" ])
    if ctx.env.ESGLOBAL == 64:
        ctx.env.append_unique("CXXFLAGS", [ "-Desglobal=long", "-Desglobal_mpi=MPI_LONG" ])

def set_default(ctx):
    try:
        if ctx.options.mpich:
            ctx.find_program("mpic++.mpich", var="MPICXX", mandatory=False)
            if not ctx.env.MPICXX:
                ctx.find_program("mpic++", var="MPICXX")
        else:
            ctx.find_program("mpic++", var="MPICXX")

    except ctx.errors.ConfigurationError:
        ctx.fatal("mpic++ not found.\nRun './waf --help' for more options.")

    if ctx.options.mpich:
        ctx.env.append_unique("CXXFLAGS", [ "-cxx=icpc" ])
        ctx.env.append_unique("LINKFLAGS", [ "-cxx=icpc" ])

    ctx.env.CXX = ctx.env.LINK_CXX = ctx.env.MPICXX;
    ctx.options.gfortran = True

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


