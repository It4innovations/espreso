
===================================
Installation of the library ESPRESO
===================================


Download the library
--------------------

ESPRESO is available through `Git <https://git-scm.com/>`_. The main repository is located at `GitLab <https://code.it4i.cz/mec059/espreso>`_.
The source codes can be easily downloaded by cloning the following repository: ::

  $ git clone git@code.it4i.cz:mec059/espreso.git
    
If you do not have an account at code.it4i.cz and do not participate in the development of the ESPRESO please use the https interface: ::

  $ git clone https:/code.it4i.cz/mec059/espreso.git

The repository contains several branches. The stable version of the library can be found in the ``master`` branch which is also the default branch when ``git clone`` is executed.

Directory structure
^^^^^^^^^^^^^^^^^^^

The main directory of the ESPRESO contains the following folders:

 - **doc** - the source code of this documentation [`ESPRESO <index.html>`__]
 - **env** - scripts related to environment settings [`Set up the environment`_]
 - **examples** - simple configuration examples for showing how to setup the ESPRESO [`examples <examples.html>`__]
 - **libespreso** - API for interfacing with 3D party applications [`ESPRESO API <api.html>`__]
 - **machines** - scripts for launching the ESPRESO on selected supercomputers [`Support for HPC <api.html>`__]
 - **src** - source files of the library [`Structure of the library <structure.html>`__]
 - **tests** - scripts for ESPRESO installation validation [`Testing the installation`_]
 - **tools** - third party tools required by the ESPRESO [`Dependencies`_]

Besides the above directories, there are also configuration and build scripts which are described in the `Building the ESPRESO`_ section.


Dependencies
------------

In the current version the ESPRESO can be compiled and executed on a Linux operating system only.
As of now, the library requires the `Intel MKL <https://software.intel.com/en-us/intel-mkl>`_ and `METIS <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`_ libraries.
Prefered and tested compiler is the `Intel compiler <https://software.intel.com/en-us/intel-compilers>`_ - another compilers are also supported, but not fully tested.

ESPRESO has interface for several Sparse Direct solvers. Depending on the selected solver (`Support for various direct solvers`_) ESPRESO requires the following libraries:

 - Pardiso (both the original version and the MKL versions are suported)  
 - MUMPS - not thread-safe solver, works only if ESPRESO is compiled without threading support
 - CuSolver - GPU accelerated solver (experimental support - low performance) 

Another mandatory system libraries:

 - pthread
 - libz
 - cmake - minimal version 2.8


Building the ESPRESO
--------------------

To compile and install the ESPRESO a Python-based framework `Waf <https://waf.io/book/>`_ is used. 
The compilation process has two phases: ``configuration`` and ``installation``. 
The former configures the persistent data and checks all required headers and libraries. 
The latter builds the library and compile the main executable file ``espreso``. 

Configuration
^^^^^^^^^^^^^

The **configuration** process is started by the following command: ::

  $ ./waf configure

It configures the library based on the default settings described in the ``build.config.default`` file. Machine specific configurations are also supported.
For instance, a configuration for ``Cray`` machines is done using: ::

  $ ./waf configure --cray

In this case, the default configuration values are taken from the ``build.config.cray`` file.

When something goes wrong and configuration crashes, you may want to check and potentialy change some of the settings. This should not be done
directly in the ``build.config.default`` (or ``build.config.cray``) file. Instead, it is recommended to create a new file ``build.config``
from the default configuration file: ::

  $ cp build.config.default build.config

When you create the ``build.config`` file, all setttings defined in this file redefine the settings from the default configuration file.
The configuration process accepts the following attributes:


 - CXX - C++ compiler (if MPI is used user needs to specifie the MPI/C++ compiler) 
 - CC - C compiler
 - FC - Fortran compiler

 + INT_WIDTH - integer width (32 - default or 64 bits - for large problems of size over 2.1 billion unknowns)
 + SOLVER - Sparse Direct Solver package (`Support for various direct solvers`_)
 + LIBTYPE - type of produced library
 + BUILD_TOOLS - ESPRESO compiles third party tools and libraries
 + METISLIB - the file name of the METIS library, if BUILD_TOOLS is set to 0

 - CXXFLAGS - general C++ compiler flags
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
The suported sparse direct solvers are: 

  - MKL
  - PARDISO
  - MUMPS


Hardware Acceleration
^^^^^^^^^^^^^^^^^^^^^
By default the CPU is used for processing. However, ESPRESO supports also modern hardware in form of accelerators. The Intel Xeon Phi (MIC) and Nvidia GPU (GPU) accelerators are suported. The options are: 

 - CPU
 - GPU
 - MIC


Installation
^^^^^^^^^^^^

After setting the sparse direct solver, hardware accelerator and successful configuration, ESPRESO can be installed by calling the following command: ::

  $ ./waf install

This command builds all source files and creates the ``espreso`` executable file.
Depending on the ``LIBTYPE``, the ``libespreso/feti4i.so`` or ``libespreso/feti4i.a``
libraries are also created during the instalation. 


Set up the environment
----------------------

Before `running <run.html>`__ the ``espreso``, following environment variables needs to be set: 

 - MKL_NUM_THREADS - in the current version it should be set to 1
 - OMP_NUM_THREADS - in the current version it should be set to 1

The last three variables should be set according to the number of CPU cores per compute node (nCores) and number of MPI processes processed per node (PPN):

 - SOLVER_NUM_THREADS - should be set to nCores/PPN
 - PAR_NUM_THREADS - should be set to nCores/PPN
 - CILK_NWORKERS - should be set to nCores/PPN

Sample environment setting files can be found in the ``env`` directory.

Testing the installation
------------------------

The installation of ESPRESO can be validated by the included set of tests which can be executed as follows: ::

  $ python tests/espreso.py

If all tests pass, ESPRESO is ready to use. Congratulations !! 
