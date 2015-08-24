
import commands
import sys
import os

sys.path.append(os.path.abspath("./python"))
from wafutils import *


#..............................................................................#
#     VERSION -> When you change the configuration, increment the value

VERSION = 10

#..............................................................................#


def configure(ctx):
    try:
        ctx.load("icpc")
    except ctx.errors.ConfigurationError:
        ctx.fatal("Install Intel compiler or load the appropriate module.")

    configure_indices_width(ctx)
    ctx.define("version", VERSION)
    ctx.write_config_header("config.h")

################################################################################
################################################################################
#                  Set MPI compiler and linker for clusters

    if ctx.options.anselm:
        ctx.env.CXX = ctx.env.LINK_CXX = ["mpic++"]
    elif ctx.options.salomon:
        ctx.env.CXX = ctx.env.LINK_CXX = ["mpiicpc"]
    elif ctx.options.titan:
        ctx.env.CXX = ctx.env.LINK_CXX = ["CC"]
    else:
        set_default(ctx)

#                   Global flags used for all libraries
#..............................................................................#

    ctx.env.append_unique("CXXFLAGS", [ "-Wall", "-openmp", "-std=c++11", "-O2" ])

################################################################################
################################################################################

    set_indices_width(ctx)

    ctx.recurse("bem")
    ctx.recurse("mesh")
    ctx.recurse("permoncube")
    ctx.recurse("solver")
    ctx.recurse("app")


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                          Options for ESPRESO

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

    opt.add_option("--cuda",
        action="store_true",
        default=False,
        help="Create application with CUDA support.")

    opt.add_option("--salomon",
        action="store_true",
        default=False,
        help="Create application for Salomon.")

    opt.add_option("--titan",
        action="store_true",
        default=False,
        help="Create application for Titan.")

    opt.add_option("--pardiso",
        action="store_true",
        default=False,
        help="Solver use pardiso library.")

    opt.add_option("--mumps",
        action="store_true",
        default=False,
        help="Solver use mumps library.")

    opt.add_option("--gfortran",
        action="store_true",
        default=False,
        help="MUMPS use libraries builder by gfortran.")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


def build(ctx):

    if not check_configuration(VERSION):
        ctx.fatal("Settings of ESPRESO have changed. Run configure first.")

    ctx(
        export_includes = "include",
        name            = "espreso_includes"
    )
    ctx.ROOT = ctx.path.abspath()

    ctx.recurse("bem")
    ctx.recurse("mesh")
    if ctx.options.mesh:
        return

    ctx.recurse("permoncube")
    if ctx.options.permoncube:
        return

    ctx.recurse("solver")
    ctx.recurse("app")





