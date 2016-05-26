
===================================
Installation of the library ESPRESO
===================================


Download the library
--------------------

ESPRESO is versioned in `Git <https://git-scm.com/>`_. The main repository is at `GitLab <https://code.it4i.cz/mec059/espreso>`_.
The source codes can be easily obtained by clone this repository: ::

  $ git clone git@code.it4i.cz:mec059/espreso.git

The repository contains several developing braneches. The last stabil version is in branch ``master``.
It is the default branch after the clone.

Directory structure
^^^^^^^^^^^^^^^^^^^

The main directory of ESPRESO contains the following directories:

 - **doc** - the source of this documentation [`ESPRESO <intex.html>`__]
 - **env** - scripts for environment settings [`Set up the environment`_]
 - **examples** - simple examples runable with ESPRESO [`examples <examples.html>`__]
 - **libespreso** - API for usage in other softwares [`ESPRESO API <api.html>`__]
 - **machines** - scripts for running ESPRESO on selected cupercomputers [`Support for HPC <api.html>`__]
 - **src** - source files of the library [`Structure of the library <structure.html>`__]
 - **tests** - scripts for testing the installation [`Testing the installation`_]
 - **tools** - third party software used by ESPRESO [`Dependencies`_]

Besides above directories, there are also building scripts described in `Building the ESPRESO`_.


Dependencies
------------

In the current version ESPRESO can be compiler only in Linux.
The main functionality is build upon the `Intel MKL <https://software.intel.com/en-us/intel-mkl>`_ and `METIS <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`_ libraries.
Prefered compiler is `Intel compiler <https://software.intel.com/en-us/intel-compilers>`_ - another compilers are also supported, but not tested.

Depends on a selected direct solver (`Support for various direct solvers`_) ESPRESO needs the following libraries:

 - Pardiso
 - MUMPS
 - Cublas, Cudard

Another mandatory system libraries:

 - pthread
 - libz
 - cmake - minimal version 2.8


Building the ESPRESO
--------------------

For compilation and instalation ESPRESO uses a Python-based framework `Waf <https://waf.io/book/>`_.
The compilation includes two phases: ``configuration`` and ``installation``.
The former configures persistend data and checks all required headers and libraries for your installation.
The latter makes the library and produces the runable binary file ``espreso``.

Configuration
^^^^^^^^^^^^^

The **configuration** is made by the command: ::

  $ ./waf configure

It configures the library based on the default settings from the file ``build.config.default``.
Configuration on a ``cray`` machine is made by the command: ::

  $ ./waf configure --cray

In this case, the configuration default values are get from the file ``build.config.cray``.

When something goes wrong, you may want to change some attributes. It should not be done
directly in the file ``build.config.default`` (or ``build.config.cray``). It is recommended to create new file ``build.config``
from the default settings: ::

  $ cp build.config.default build.config

When you create the file ``build.config``, all setttings from this file rewrite attributes
from the default file.
The configuration process accepts the following attributes:


 - CXX - C++/MPI compiler:
 - CC - C compiler
 - FC - Fortran compiler

 + INT_WIDTH - integer width
 + SOLVER - direct solver (`Support for various direct solvers`_)
 + LIBTYPE - type of produced library
 + BUILD_TOOLS - ESPRESO compiles third party libraries
 + METISLIB - the name of METIS library, if BUILD_TOOLS is set to 0

 - CXXFLAGS - general compiler flags
 - LINKFLAGS - general linker flags
 - INCLUDES - general include directories
 - LIBPATH - general dynamic library paths
 - STLIBPATH - general static library paths

 + {SOLVER}::CXXFLAGS - solver specific compiler flags
 + {SOLVER}::LINKFLAGS - solver specific linker flags
 + {SOLVER}::INCLUDES - solver specific include directories
 + {SOLVER}::LIBPATH - solver specific dynamic library paths
 + {SOLVER}::STLIBPATH - solver specific static library paths


Support for various direct solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MKL, PARDISO, MUMPS


Hardware Acceleration
^^^^^^^^^^^^^^^^^^^^^
CPU, GPU, MIC


Installation
^^^^^^^^^^^^

After setting appropriate direct solver, hardware acceleration and successful configuration, ESPRESO is ready to install.
Installing is done by the command: ::

  $ ./waf install

It builds all source files and produces executable file ``espreso``.
Depend on ``LIBTYPE``, during installation is also created library ``libespreso/feti4i.so``
of ``libespreso/feti4i.a``.


Set up the environment
----------------------

Before `run <run.html>`__ the ``espreso``, the environment variables has to be set.
Sample settings files are in directory ``env``.
The following environment variables has to be set:

 - MKL_NUM_THREADS
 - OMP_NUM_THREADS
 - SOLVER_NUM_THREADS
 - PAR_NUM_THREADS

Testing the installation
------------------------

The installation can be simply tested by: ::

  $ python tests/espreso.py

If all tests pass, ESPRESO is ready to use.
