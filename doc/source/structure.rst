

Structure of the library
========================

This section describes the structure of the source files.
Source files are located in the ``src`` directory which has the following structure:

 - **app** - final applications and factory for loading the problem from input parameters
 - **assembler** - matrix assemblers for the solver
 - **basis** - general classes for parameters loading, logging, etc.
 - **catalyst** - Paraview Catalyst interface
 - **config** - the main ESPRESO configuration file
 - **include** - headers of the third party libraries
 - **input** - ESPRESO mesh loaders
 - **mesh** - main classes for mesh processing and problem description
 - **output** - classes for saving the mesh from ESPRESO
 - **python** - python scripts
 - **solver** - ESPRESO FETI solver

Except for the **solver**, all directories contains general code that is used by all configurations.
The requested type of solver is built based on the `build configuration <installation.html#configuration>`_ which is defined by the user.
ESPRESO adjusts everything else. The next section describes how this feature is implemented.

Implementation of solver switching
__________________________________

The implementation can be divided into two parts:
``set up the environment`` and ``compilation of selected classes``.

Set up the environment
^^^^^^^^^^^^^^^^^^^^^^

The environment is set by ``waf``.
It checks the availability of headers and libraries for the selected solver.
The check is done by the following set of commands in the ``wscript``: ::

  def configure(ctx):
      ctx.check_header("header")
      ctx.check_stlib("static_lib")
      ctx.check_lib("shared_lib")


Compilation of selected classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Classes for the selected solver are defined in the ``src/solver/wscript`` file.
The ``wscript`` file contains the configuration method with following rules: ::

    if ctx.env.SOLVER == "SOLVER_TYPE":
        ctx.env.append_unique("LIB", [ "library1" ])
        ctx.env.append_unique("POSTSTLIB", [ "stlibrary3" ])
        ctx.env.append_unique("STLIB", [ "stlibrary1", "stlibrary2" ])

Where ``SOLVER_TYPE`` is value of parameter ``SOLVER`` from `build configuration <installation.html#configuration>`_,
``LIB``, ``POSTSTLIB``, and ``STLIB`` are libraries from the following groups:

  :STLIB: static libraries
  :LIB: shared libraries
  :POSTSTLIB: static libraries linked last

The compilation command generated from the above examples is: ::

  $ mpic++ FLAGS SOURCES -Wl,-Bstatic -Wl,--start-group stlibrary1 stlibrary2 -Wl,--end-group -Wl-Bdynamic library1 -Wl,-Bstatic -Wl,--start-group stlibrary3 -Wl,--end-group -Wl-Bdynamic

.. note::
   The order of ``LIB``, ``STLIB``, and ``POSTSTLIB`` is always the same (``STLIB`` - ``LIB`` - ``POSTSTLIB``).

   The libraries should also be adjusted for both static and dynamic ``LIBTYPE``.

The ``wscript`` file also contains the build method with similiar rules: ::

    if ctx.env.SOLVER == "SOLVER_TYPE":
        sources = source_files + ("file1", "file2", "file2")

By this rule, the appropriate source files are compiled.
The list has to contain all the files needed by the selected solver.

At this point the build framework should be correctly set.
The next step is to prepare the content of the source files.

The solver is not aware of the selected solver type.
The solver only distinguishes between the CPU and the accelerated version.
The particular type is set in the following files:

 - src/solver/specific/clusters.h
 - src/solver/specific/densesolvers.h
 - src/solver/specific/sparsesolvers.h
 - src/solver/specific/itersolvers.h

The files contain typedefs for generic types.
For example the sparsesolvers.h has the following tepedefs: ::

   #if defined(SOLVER_MKL)
   #include "cpu/SparseSolverMKL.h"

   namespace espreso {
      typedef SparseSolverMKL SparseSolverCPU;
      typedef SparseSolverMKL SparseSolverAcc;
   }

The meaning of the above example is:
if ``SOLVER==MKL`` then include ``SparseSolverMKL`` and use it for both CPU and Accelerated version.








