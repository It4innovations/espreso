
import commands
import os

VERSION = 2

def options(opt):
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

    opt.add_option("--permoncube",
       action="store_true",
       default=False,
       help="Create application only from permoncube.")

    opt.add_option("--anselm",
       action="store_true",
       default=False,
       help="Create application for Anselm.")


def configure(ctx):
    if ctx.options.anselm:
        anselm(ctx)
    else:
        check_environment(ctx)

    data_types(ctx)
    write_configuration(ctx)

    ctx.setenv("base", ctx.env)
    ctx.env.append_unique("CXXFLAGS", [ "-Wall", "-openmp" ])
    ctx.env.append_unique("LINKFLAGS", [ "-Wall", "-openmp" ])

    if ctx.env.ESINT == 32:
        ctx.env.append_unique("CXXFLAGS", [ "-Desint=int", "-DMKL_INT=int" ])
        ctx.env.append_unique("LIB", [ "mkl_intel_lp64" ])
    if ctx.env.ESINT == 64:
        ctx.env.append_unique("CXXFLAGS", [ "-Desint=long", "-DMKL_INT=long" ])
        ctx.env.append_unique("LIB", [ "mkl_intel_ilp64" ])

    if ctx.env.ESLONG == 32:
        ctx.env.append_unique("CXXFLAGS", [ "-Deslong=int" ])
    if ctx.env.ESLONG == 64:
        ctx.env.append_unique("CXXFLAGS", [ "-Deslong=long" ])

    ctx.env.append_unique("LIB", [ "mkl_core" ])

    if ctx.options.debug:
        ctx.env.append_unique("CXXFLAGS", [ "-g" ])
        ctx.env.append_unique("LIB", [ "mkl_sequential" ])
    else:
        ctx.env.append_unique("CXXFLAGS", [ "-O2", "-DXE6", "-DDEVEL", "-DTM_BLOCK_START"])
        ctx.env.append_unique("LIB", [ "mkl_intel_thread" ])

    ctx.env.append_unique("LIB", [ "pthread", "m" ])

    ctx.setenv("mpi", ctx.env)
    if ctx.options.mpich:
        ctx.env.append_unique("CXXFLAGS", [ "-cxx=icpc" ])
        ctx.env.append_unique("LINKFLAGS", [ "-cxx=icpc" ])

    if ctx.options.anselm:
        ctx.env.append_unique("CXXFLAGS", [ "-xSSE4.1" ])

    ctx.env.CXX = list(ctx.env.MPICXX)
    ctx.env.LINK_CXX = list(ctx.env.MPICXX)

    ctx.recurse("metis")
    ctx.recurse("bem")
    ctx.recurse("mesh")
    ctx.recurse("permoncube")
    ctx.recurse("solver")
    ctx.recurse("app")

def build(ctx):
    test_file = os.path.join(ctx.path.abspath(), "build/test_config")
    if not os.path.isfile(test_file):
        write_test_file()

    process, version = commands.getstatusoutput(test_file)
    if int(version) != VERSION:
        ctx.fatal("Settings of ESPRESO have changed. Run configure first.")

    ctx.ROOT = ctx.path.abspath()

    ctx.recurse("metis")
    ctx.recurse("bem")
    ctx.recurse("mesh")
    if ctx.options.mesh:
        return

    ctx.recurse("permoncube")
    if ctx.options.permoncube:
        return

    ctx.recurse("solver")
    ctx.recurse("app")



def anselm(ctx):
    ctx.load("icpc")
    ctx.env.MPICXX = ["mpic++"]


def check_environment(ctx):
    try:
        ctx.load("icpc")
    except ctx.errors.ConfigurationError:
        ctx.fatal("Install Intel compiler or try configuration for your cluster.\n"
            "Run './waf --help' for more options.")

    ctx.check(
        "cmake",
        msg = "Checking for cmake"
    )

    try:
        if ctx.options.mpich:
            ctx.find_program("mpic++.mpich", var="MPICXX")
        else:
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
        define_name = "esint",
        errmsg      = "Incorrect user-supplied value for ESPRESO_LOCAL_INDICES_WIDTH",
        msg         = "Checking ESPRESO_LOCAL_INDICES_WIDTH"
    )
    esint = ctx.get_define("esint").replace("\"", "")
    ctx.define("esint", int(esint))
    ctx.env.ESINT = int(esint)

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
        define_name = "eslong",
        errmsg      = "Incorrect user-supplied value for ESPRESO_GLOBAL_INDICES_WIDTH",
        msg         = "Checking ESPRESO_GLOBAL_INDICES_WIDTH"
    )
    eslong = ctx.get_define("eslong").replace("\"", "")
    ctx.define("eslong", int(eslong))
    ctx.env.ESLONG = int(eslong)


def write_configuration(ctx):
    ctx.define("version", VERSION)
    ctx.write_config_header("config.h")
    write_test_file()


def write_test_file():
    test_config = open("build/test_config", "w")
    test_config.write('#!/usr/bin/env python')
    test_config.write(
    '''
import os

if not os.path.isfile("build/config.h") or not os.path.isfile("include/espreso.h"):
    print "0"
    quit()

configFile = open("build/config.h", "r")
espresoFile = open("include/espreso.h", "r")

config = {}
espreso = {}

for line in configFile:
    if line.startswith("#define"):
        values = line.split()
        config[values[1]] = values[-1]
for line in espresoFile:
    if line.startswith("#define"):
        values = line.split()
        espreso[values[1]] = values[-1]
if config["esint"] != espreso["ESPRESO_LOCAL_INDICES_WIDTH"]:
    print "0"
elif config["eslong"] != espreso["ESPRESO_GLOBAL_INDICES_WIDTH"]:
    print "0"
else:
    print config["version"]
    ''')

    os.system("chmod +x build/test_config")
