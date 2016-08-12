

Run the ESPRESO
===============

After the successful `installation <installation.html>`__
and setting the `environment <installation.html#set-up-the-environment>`__,
there is an executable file ``espreso`` in the main directory.
The following sections describe how to use and execute the ESPRESO.

Settings input parameters
-------------------------

The execution of the ESPRESO consists of two phases: 

  - generating the `problem <run.html#specification-of-a-problem>`_ based on the problem description file, 
  - running the FETI `solver <run.html#proper-setting-of-solver>`_ to calculate the solution.

Both phases are controled by the input parameters.
The `espreso` file takes input parameters in two ways: (i)
from configuration file, and (ii) from command line arguments.
The basic list of input parameters can be listed by calling: ::

  $ ./espreso -h

The basic parameters, listed by the command above, should be sufficient for standard usage of the ESPRESO.
The experimental settings for advanced users and develpers can be found in the exhaustive list of parameters: ::

  $ ./espreso -hh

Please note: Parameters from this list adjust internal settings that are not fully tested
or are under development. Therefore, setting the wrong values can lead to the unexpected behaviour.

The prefered way to set the input parameters of the ESPRESO is using the **Configuration file**.
By default, ``espreso`` uses ``espreso.config`` configuration file located in the same directory as the executable file.
Different configuration files can be specified as follows: ::

  $ ./espreso -c arbitrary_config_file

The ESPRESO configuration file contains a list of parameters and assigned values in the following form:

``PARAMETER = VALUE``

where ``PARAMETER`` is the name of an input parameters (`complete list <parameters.html>`__).
The ``PARAMETER`` name is case-insensitive. 
Configuration file uses ``#`` characters for comments, therefore everything following this character in the current line is ignored.
The ``VALUE`` can be defined using one of the following types:

 :INTEGER:  *integer value*
 :DOUBLE:   *floating point value*
 :STRING:   *case-sensitive string value*
 :BOOLEAN:  *value* ``0`` *is interpreted as* ``false`` *, everything else (including empty value) is interpreted as* ``true``

**Command line arguments** provide alternative way of setting the input parameters.
They override the settings from the configuration file.
The following arguments are accepted:

 -h                  basic help message (-hh for advanced help)
 -i TYPE             sets the type of problem description/input data: TYPE={workbench, openfoam, generator, esdata}
 -p PATH             sets a path to an problem description file
 -c FILE             sets a file with the ESPRESO settings
 -v                  increases the verbose level (can be applied recursivelly, e.g. -vv)
 -t                  increases the testing level
 -m                  increases the processing time measurement verbose level
 --PARAMETER=VALUE   overrides a value of one particular parameter from the configuration file (*Note: the name is case-sensitive!*)

Specification of a problem
--------------------------

ESPRESO supports a problem description in two widely used formats: `Ansys Workbench <http://www.ansys.com/>`__ and `OpenFOAM <http://www.openfoam.com/>`__.
In addition, ESPRESO have a general `API <api.html>`__ that allows third party applications to use ESPRESO as a parallel linear solver.
For testing and edvelopment purposes ESPRESO also contains a problem generator.


Ansys
^^^^^

ESPRESO can open database files in `Ansys Workbench <http://www.ansys.com/Products/Platform>`__ format and solve the problems they describe. If the problem is relativelly small, round 5 million unknowns, it can be solved on a single compute node/workstation with 64 GB of RAM. 

ESPRESO is designed to solve large problems using supercomputers with many compute nodes. To run the solver on multiple compute nodes the original database file needs to decomposed (using the domain decomposition approach) into multiple database files. 

For this decomposition, the ESPRESO contains a ``decomposer`` tool. The ``depomposer`` is automatically compiled during the installation process.
It is executed on a single compute node or workstation preferably with large amount of main memory. It loads the Ansys database file and save the decomposed problem into files using the ESPRESO internal binary format suitable for parallel loading.

 
The ``decomposer`` accepts has following input parameters:

  - a path including the file name to the input batabase file (here: ./workbench.dat)
  - a path to an output directory (here: ./decomposition)
  - a number of output files (should be the same as number of MPI processes). Please note: one can specify multiple decompositions that will be carried out at the same time, see the next example, which decomposes the input problem into 32, 64 and 128 subdomains. This way one can prepare data for scalability tests.

::

  $ ./decomposer workbench.dat decomposition 32 64 128

This example loads ``workbench.dat`` and creates three directories: ``decomposition32``, ``decomposition64``, and ``decomposition128``.
The directories contain the ``workbench.dat`` example decomposed into 32, 64, and 128 parts. 

The ``decopomser`` can also be used for a hierarchical decomposition of a previously decomposed problem. The parameters are the same as in the previous example.
However, the ``decomposer`` needs to be executed in parallel: ::

  $ mpirun -n 128 ./decomposer decomposition128 decomposition 4 8

The above example creates directories ``decomposition512`` (128 * 4) and ``decomposition1024`` (128 * 8).

.. note::
   The parallel decomposition does not re-construct the original problem.
   The decomposition is executed using the previously decomposed parts.
   Hence, multi-level decomposition is **not the same** as the single-level decomposition,
   even though the final number of parts is the same!
   It is strongly recommended to use single decomposition if possible
   to minimize the load balancing problems.

Finally the decomposed problem can be solved by following command: ::

  $ mpirun -n 512 ./espreso -i esdata -p decomposition512

Where ``-i esdata`` specifies the input data format (esdata is a format used by the decomposer) and ``-p decomposition512`` defines the input directory. 

OpenFOAM
^^^^^^^^

The OpenFOAM support is under development.

Problem generator
^^^^^^^^^^^^^^^^^

Is a tool which generates a mesh with an arbitrary number of elements of particular type.
It is the key tool for the solver scalability tests on the masivelly parallel machines as it is able to generate a multi-billion problem in several seconds. 
Internally, it is also widely used for the development and testing of the new features.

Scripts with generator settings can be found in ``examples/meshgenerator/`` directory.
Scripts are in the same format as the ESPRESO configuration files (``PARAMETER = VALUE``)
and usualy accepts list of nameless command line parameters.
These parameters usualy specify the size of the generated problem and the element type.
For example, the following command: ::

  $ mpirun -n 8 ./espreso -p examples/meshgenerator/cube_elasticity_fixed_bottom.txt HEXA8  4 2 1  2 4 8  2 2 4

generates and solves a cubical problem composed of:

  - hexahedron elements (HEX8), 
  - the problem is divided into 4*2*1 = 8 MPI processes (4 in X direction, 2 in Y direction and 1 in Z direction),
  - problem on each MPI process is further decomposed into to 2*4*8 = 64 subdomains (2 in X direction, 4 in Y direction and 8 in Z direction)
  - each subdomain is composed of 2*2*4 = 16 elements (2 in X direction, 2 in Y direction and 4 in Z direction)


Check the results
-----------------

The solution can be viewed in `Paraview <http://www.paraview.org/>`__.
By default the ``espreso`` save results into legacy VTK files using the following naming convention ``result{MPI_rank}.vtk``.



