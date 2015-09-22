
import commands
import sys
import os

sys.path.append(os.path.abspath("./python"))
from wafutils import *


#..............................................................................#
#     VERSION -> When you change the configuration, increment the value

VERSION = 13

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
        ctx.env.CXX = ctx.env.LINK_CXX = [ "mpic++" ]
        ctx.env.MPICC = [ "mpiicc" ]
        ctx.env.MPIFORT = [ "mpiifort" ]
    elif ctx.options.salomon:
        ctx.env.CXX = ctx.env.LINK_CXX = ["mpiicpc"]
        ctx.env.MPICC = [ "mpiicc" ]
        ctx.env.MPIFORT = [ "mpiifort" ]
    elif ctx.options.titan:
        ctx.env.CXX = ctx.env.LINK_CXX = ["CC"]
        ctx.env.MPICC = [ "cc" ]
        ctx.env.MPIFORT = [ "fc" ]
    else:
        set_default(ctx)

#                   Global flags used for all libraries
#..............................................................................#

    ctx.env.append_unique("CXXFLAGS", [ "-Wall", "-openmp", "-std=c++11", "-O2" ])
    ctx.env.append_unique("LINKFLAGS", [ "-Wall", "-openmp" ])
    if ctx.options.titan:
        ctx.env.append_unique("CXXFLAGS", [ "-fPIE", "-dynamic" ])
        ctx.env.append_unique("LINKFLAGS", [ "-pie", "-dynamic" ])


################################################################################
################################################################################

    set_indices_width(ctx)

    ctx.ROOT = ctx.path.abspath()

    ctx.recurse("tools")
    ctx.recurse("bem")
    ctx.recurse("mesh")
    ctx.recurse("permoncube")
    ctx.recurse("solver")
    ctx.recurse("catalyst")
    ctx.recurse("app")


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                          Options for ESPRESO


def options(opt):

# Configure options

    opt.add_option("--anselm",
        action="store_true",
        default=False,
        help="Create application for Anselm.")

    opt.add_option("--salomon",
        action="store_true",
        default=False,
        help="Create application for Salomon.")

    opt.add_option("--titan",
        action="store_true",
        default=False,
        help="Create application for Titan.")

    opt.add_option("--mpich",
        action="store_true",
        default=False,
        help="Compile MPI version with mpich.")

    opt.add_option("--gfortran",
        action="store_true",
        default=False,
        help="MUMPS use libraries builder by gfortran.")


# Build options

    opt.add_option("--static",
        action="store_true",
        default=False,
        help="All libraries created by ESPRESO are static.")

    opt.add_option("--mesh",
        action="store_true",
        default=False,
        help="Create application only from mesh.")

    opt.add_option("--cuda",
        action="store_true",
        default=False,
        help="Create application with CUDA support.")

    opt.add_option("--mic",
        action="store_true",
        default=False,
        help="Create application with MIC support.")

    opt.add_option("--pardiso",
        action="store_true",
        default=False,
        help="Solver use pardiso library.")

    opt.add_option("--mumps",
        action="store_true",
        default=False,
        help="Solver use mumps library.")

    opt.add_option("--catalyst",
        action="store_true",
        default=False,
        help="ESPRESO supports ParaView-Catalyst.")
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


def build(ctx):

    if not check_configuration(VERSION):
        ctx.fatal("Settings of ESPRESO have changed. Run configure first.")

    ctx(
        export_includes = "include",
        name            = "espreso_includes"
    )
    ctx.ROOT = ctx.path.abspath()
    ctx.LIBRARIES = ctx.ROOT + "/libs"

    if ctx.options.static:
        ctx.lib = ctx.stlib
    else:
        ctx.lib = ctx.shlib

    if ctx.options.mic:
        pass
        #ctx.recurse("solver")
        #return

    ctx.recurse("tools")
    ctx.add_group()
    ctx.recurse("bem")
    ctx.recurse("mesh")
    if ctx.options.mesh:
        return

    ctx.recurse("permoncube")
    ctx.recurse("solver")
    if ctx.options.catalyst:
        ctx.recurse("catalyst")
    ctx.recurse("app")





