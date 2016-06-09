

Structure of the library
========================

This section describes the structure of source files.
Source files are in directory ``src`` and have the following structure:

 - **app** - final applications and factory for loading the problem from input parameters
 - **assembler** - assemblers of matrices for the solver
 - **basis** - general classes for load parameters, logging, etc.
 - **catalyst** - Paraview Catalyst wrapper
 - **config** - the main ESPRESO configuration file
 - **include** - headers of third parties libraries
 - **input** - ESPRESO mesh loaders
 - **mesh** - main classes for description of a problem
 - **output** - classes for saving the ESPRESO mesh
 - **python** - python scripts
 - **solver** - ESPRESO solver

Except **solver**, all directories contains general code that is used in all configurations.
An appropriate solver is builded based on the `build configuration <installation.html#configuration>`_.
A user is responsible only for selection a correct solver on a given architecture.
ESPRESO adjust everything else. The next section describes how this feature is implemented.

Implementation of compilation an appropriate solver
___________________________________________________

The implementation can be divided into the two parts:
``set up the environment`` and ``compilation of selected classes``.

Set up the environment
^^^^^^^^^^^^^^^^^^^^^^

The environment is set by ``waf``.
It checks the availability of headers and libraries for chosen solver.
Checking is based on the settings in source file ``src/python/wafutils.py``.
The file contains method ``check_libraries``, where the libraries are listed.

.. warning::
   TODO: check will be refactored


Compilation of selected classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Classes for a selected solver are set in file ``src/solver/wscript``.
The ``wscript`` file contains the configuration method with rules of type: ::

    if ctx.env.SOLVER == "SOLVER_TYPE":
        ctx.env.append_unique("LIB", [ "library1" ])
        ctx.env.append_unique("POSTSTLIB", [ "stlibrary3" ])
        ctx.env.append_unique("STLIB", [ "stlibrary1", "stlibrary2" ])

Where ``SOLVER_TYPE`` is value of parameter ``SOLVER`` from `build configuration <installation.html#configuration>`_,
``LIB``, ``POSTSTLIB``, and ``STLIB`` are constants separated libraries into three groups:

  :STLIB: static libraries
  :LIB: shared libraries
  :POSTSTLIB: static libraries linked last

The compilation command generated from the above examples is: ::

  $ mpic++ FLAGS SOURCES -Wl,-Bstatic -Wl,--start-group stlibrary1 stlibrary2 -Wl,--end-group -Wl-Bdynamic library1 -Wl,-Bstatic -Wl,--start-group stlibrary3 -Wl,--end-group -Wl-Bdynamic

.. note::
   The order of ``LIB``, ``STLIB``, and ``POSTSTLIB`` is always the same (``STLIB`` - ``LIB`` - ``POSTSTLIB``).

   The libraries should also be adjusted for both static and dynamic ``LIBTYPE``.

The ``wscript`` file also contains the build method the simillar rules: ::

    if ctx.env.SOLVER == "SOLVER_TYPE":
        sources = source_files + ("file1", "file2", "file2")

By this rule, appropriate source files are added to sources for compilation.
The list has to contains all files needed by the selected solver.

At this point the build framework should be correctly set.
The next step prepares the source files.

The solver is not aware of the selected solver type.
The solver only distinguishes between CPU or accelerated version,
and a particular type is set in the following files:

 - src/solver/specific/clusters.h
 - src/solver/specific/densesolvers.h
 - src/solver/specific/sparsesolvers.h
 - src/solver/specific/itersolvers.h

The files contains typedefs for generic types.
For example the sparsesolvers.h has the following tepedefs: ::

   #if defined(SOLVER_MKL)
   #include "cpu/SparseSolverMKL.h"

   namespace espreso {
      typedef SparseSolverMKL SparseSolverCPU;
      typedef SparseSolverMKL SparseSolverAcc;
   }

The meaning of the above example is:
if ``SOLVER==MKL`` then include ``SparseSolverMKL`` and use it for both CPU and Accelerated computation.








