
===================================
Installation of the library ESPRESO
===================================


Download the library
--------------------

ESPRESO is available through `Git <https://git-scm.com/>`_. The release repository is located at `GitLab <https://github.com/It4innovations/espreso>`_.
The source codes can be downloaded by cloning the following repository (for users with access to private IT4I repositories): ::

  $ git clone git@code.it4i.cz:mec059/espreso.git
    
If you do not have an account at code.it4i.cz and do not participate in the development of the ESPRESO please use the https interface: ::

  $ git clone https://github.com/It4innovations/espreso

The repository contains several branches. The stable version of the library can be found in the ``master`` branch which is also the default branch when ``git clone`` is executed.

Directory structure
^^^^^^^^^^^^^^^^^^^

The main directory of the ESPRESO contains the following folders:

 - **benchmarks** - simple benchmarks for showing how to setup the ESPRESO [`examples <examples.html>`__]
 - **doc** - the source code of this documentation [`ESPRESO <index.html>`__]
 - **env** - scripts related to environment settings [`Set up the environment`_]
 - **install** - configuration files for building the library [`Building the ESPRESO`_]
 - **libespreso** - API for interfacing with 3rd party applications [`ESPRESO API <api.html>`__]
 - **machines** - scripts for launching the ESPRESO on selected supercomputers [`Support for HPC <hpc.html>`__]
 - **src** - source files of the library [`Structure of the library <structure.html>`__]
 - **tests** - scripts for ESPRESO installation validation [`Testing the installation`_]
 - **tools** - third party tools required by the ESPRESO [`Dependencies`_]


Dependencies
------------

In the current version the ESPRESO can be compiled and executed on a Linux operating system only.
As of now, the library requires the `Intel MKL <https://software.intel.com/en-us/intel-mkl>`_ and `METIS <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`_ libraries.
Intel MKL should be installed in the system. METIS is inscluded in ESPRESO reposity and is builded automatically.
Prefered and tested compiler is the `Intel compiler <https://software.intel.com/en-us/intel-compilers>`_ or GCC.

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
If something in the environment is changed, it is mandatory to re-run the configuration process.

Configuration
^^^^^^^^^^^^^

The configuration process is controlled by ``build.config`` script.
The prefered way, how to create this script is to copy it from ``install`` directory.
The following command sets ESPRESO for compiling with Intel: ::

$ cp install/build.config.icpc build.config

The script contains the list of attributes that influence the installation: ::

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


When all attributes are correctly set, the compilation should be started by the command: ::

  $ ./waf configure


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
 - CILK_NWORKERS - in the current version it should be set to 1


The last three variables should be set according to the number of CPU cores per compute node (nCores) and number of MPI processes processed per node (PPN):

 - OMP_NUM_THREADS - should be set to nCores/PPN
 - SOLVER_NUM_THREADS - should be set to nCores/PPN
 - PAR_NUM_THREADS - should be set to nCores/PPN


Environment setting files can be found in the ``env`` directory.
The simples way how to sets the environemnt is: ::

  $ . env/paths.defaul  (sets LD_LIBRARY_PATH and PATH)
  $ . env/threading.defaul nCores/PPN  (sets above environement variables)

Testing the installation
------------------------

The installation of ESPRESO can be validated by the included set of benchmarks which can be executed as follows: ::

  $ python tests/benchmarks.py

If all tests pass, ESPRESO is ready to use. Congratulations !! 

