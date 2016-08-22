

Examples
========

In this section a set of procedures for instalation and execution of the solver is presented. 

Installation of the ESPRESO
---------------------------

Installation procedure for selected machine configurations.

CPU version with 32-bit integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the default settings for ESPRESO and no changes of the build scripts are required. ::

  for users: $ git clone https://code.it4i.cz/mec059/espreso.git
  for developers: $ git clone git@code.it4i.cz:mec059/espreso.git
  $ cd espreso
  $ ./waf configure
  $ ./waf install

CPU version with 64-bit integers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The first part is the same as for the 32-bit version::

  for users: $ git clone https://code.it4i.cz/mec059/espreso.git
  for developers: $ git clone git@code.it4i.cz:mec059/espreso.git
  $ cd espreso

ESPRESO is a 32-bit solver by default.
To enable the 64-bit integer support, the ``build.config`` has to be modified (See `configuration <installation.html#configuration>`__). ::

  $ cp build.config.default build.config

Now change the INT_WIDTH value in the build.config to: INT_WIDTH = 64, configure and install the ESPRESO: ::

  $ ./waf configure
  $ ./waf install


Enabling the Intel Xeon Phi Acceleration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This procedure is designed for the IT4Innovations Salomon cluster.  

Enabling the hardware acceleration of the ESPRESO on a cluster with accelerators requires these steps:
  1. Load the appropriate modules - set the environment 
  2. Change the build.config script
  3. Configure and install

The following example shows the installation with ``Intex Xeon Phi x100 models`` support on the ``IT4Innovations Salomon`` cluster. 
The user has to connect to the accelerated node in order to be able to compile the Xeon Phi support (using the PBS and the qsub tool). On the compute node the ESPRESO can be build as follows: ::

  $ cd espreso (or whereever the ESPRESO is cloned from the Git)
  $ . env/modules.salomon
  $ cp build.config.default build.config

Now change the value of the ``SOLVER`` in the build.config to: SOLVER = MIC and configure and build the ESPRESO: ::

  $ ./waf configure
  $ ./waf install

.. warning ::
  ``env/modules.salomon`` can be out of date.
  It is recommended to check the modules on Salomon and change the modules in this script to the latest version

.. warning ::
  ``Salomon`` has incorrectly set MIC_LD_LIBRARY_PATH.
  The correct path is: export MIC_LD_LIBRARY_PATH=/apps/compiler/icc/2016.1.150-GCC-4.9.3/mkl/lib/mic/:/apps/compiler/icc/2016.1.150-GCC-4.9.3/lib/mic/

Installing on Cray
^^^^^^^^^^^^^^^^^^

The first step is loading the required modules to setup the environment.
The following example shows the installation on ``CSC Sisu`` Cray XC40 supercomputer: ::

  $ . env/modules.sisu

  for users: $ git clone https://code.it4i.cz/mec059/espreso.git
  for developers: $ git clone git@code.it4i.cz:mec059/espreso.git
  $ cd espreso
  $ ./waf configure --cray
  $ ./waf install


Using ESPRESO as an Elmer solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following examples assumes that ESPRESO is already installed: ::

  $ mkdir elmer
  $ cd elmer
  $ git clone https://code.it4i.cz/mec059/elmer.git src

  $ mkdir build
  $ FETI4I_ROOT=${PATH_TO_ESPRESO}/libespreso cmake -DWITH_ELMERGUI:BOOL=FALSE -DWITH_MPI:BOOL=TRUE -DWITH_FETI4I:BOOL=TRUE -DCMAKE_INSTALL_PREFIX=../ ../src/
  $ make install

  $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PATH_TO_ELMER}/elmer/lib/elmersolver
  $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PATH_TO_ESPRESO}/libespreso
  $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PATH_TO_ESPRESO}/libs
  $ export PATH=$PATH:${PATH_TO_ELMER}/bin
  $ cd fem/tests/WinkelNavierPartitionUniform/
  $ ElmerGrid 1 2 winkel.grd -partition 2 2 2 2
  $ mpirun -n 4 ElmerSolver


Run the solver
--------------

All examples assume that the environment is set (See `environment settings <installation.html#set-up-the-environment>`__).


Linear elasticity problem from meshgenerator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example generates a cube that is fixed on the bottom plane and only the gravity force is applied ::

  $ mpirun -n 8 ./espreso -p examples/meshgenerator/cube_elasticity_fixed_bottom.txt HEXA8 2 2 2  5 5 5  8 8 8

Where: 
  - generator creates 8 (2x2x2) clusters
  - each cluster contains 125 (5x5x5) subdomains
  - each subdomain constains 512 (8x8x8) hexahedron elements

Detailed description of the generator parameters can be found in the ``cube_elasticity_fixed_bottom.txt`` example file.

Other examples for generator are in the ``examples/meshgenerator/`` directory.


Ansys Workbench example
^^^^^^^^^^^^^^^^^^^^^^^

The Ansys Workbench database file can be solved by one MPI process only.
To use more compute nodes the problem has to be decomposed into multiple parts. Then the solver can run in parallel: ::

  $ ./decomposer workbench_test_case.dat decomposition 4
  $ mpirun -n 4 ./espreso -i esdata -p decomposition4/







