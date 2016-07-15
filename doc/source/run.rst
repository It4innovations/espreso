

Run the ESPRESO
===============

Since the successful `installation <installation.html>`__
and setting the `environment <installation.html#set-up-the-environment>`__,
there is runable file ``espreso`` in the main directory.
The next sections constains information how to manipulate with this file.

Settings input parameters
-------------------------

The run of ESPRESO consists of two phases.
At first, the `problem <run.html#specification-of-a-problem>`_ specification is loaded, and then
the `solver <run.html#proper-setting-of-solver>`_ produces the solution.
The both phases are set by input parameters.
The `espreso` accepts input parameters in two form -
via configuration file or by command line arguments.
The basic list of input parameters can be listed by: ::

  $ ./espreso -h

Basic parameters should be sufficient for standard usage of ESPRESO.
For advanced users an exhaustive list can be listed be: ::

  $ ./espreso -hh

Parameters from this list adjust internal settings that are not fully tested
or are under development. Hence, wrong values can lead to unexpected behaviour.

**Configuration file** serves the prefered way of settings input parameters.
In the default, ``espreso`` assumes configuration file with name ``espreso.config``.
This name can be changed to arbitrary file: ::

  $ ./espreso -c arbitrary_config_file

The file consists of a number of assignments of the form:

``PARAMETER = VALUE``

where ``PARAMETER`` is one of the listed input parameters (`complete list <parameters.html>`__).
Everything after character ``#`` is omited.
The parameters are case-insensitive, values are of the following types:

 :INTEGER:  *integer value*
 :DOUBLE:   *floating point value*
 :STRING:   *case-sensitive string value*
 :BOOLEAN:  *value* ``0`` *is interpreted as* ``false`` *, everything else (including empty value) as* ``true``

**Command line arguments** provide alternative way of settings input parameters.
They override settings from the configuration file.
The following arguments are accepted:

 -h                  show basic help message (-hh for advanced help)
 -i TYPE             set type of an example TYPE={workbench, openfoam, generator, esdata}
 -p PATH             set path to an input example
 -c FILE             file with the ESPRESO configuration
 -v                  increase verbose level (can be applied more times, e.g. -vv)
 -t                  increase testing level
 -m                  increase measuring level
 --PARAMETER=VALUE   overrides value of the parameter from the configuration file (*note: name is case-sensitive!*)

Specification of a problem
--------------------------

ESPRESO supports specification of a problem in two widely used formats: `Ansys <http://www.ansys.com/>`__ and `OpenFOAM <http://www.openfoam.com/>`__.
In addition, ESPRESO provides general `API <api.html>`__ for utilize its solver in other softwares.
For testing purposes ESPRESO also contains problem generator.
It can be used for proper settings of environment etc.

Ansys
^^^^^

ESPRESO supports `Ansys Workbench <http://www.ansys.com/Products/Platform>`__ format for a problem specification.
The file should contains the geometry and boundary conditions of an example.
Because ESPRESO aims to usage on the largest supercomputers,
it is not appropriate to load a specification from a single file.
Hence, ESPRESO provides ``decomposer`` to decompose whole problem to smaller parts.
The ``depomposer`` is automatically created during the installation process.
It loads an Ansys file and save a decomposition to the ESPRESO internal binary format suitable for highly paralell loading.
The ``decomposer`` accepts a path to an input,
a path to an output,
and an array of numbers - MPI processes.
For example: ::

  $ ./decomposer workbench.dat decomposition 32 64 128

loads ``workbench.dat`` and creates three directories: ``decomposition32``, ``decomposition64``, and ``decomposition128``.
Directories contain the ``workbench.dat`` example decomposed to 32, 64, and 128 parts, respectively.
The ``decopomser`` can also be used for another decomposition of an already decomposed problem.
Parameters are the same as in the previous example.
However, the ``decomposer`` needs to be run in parallel: ::

  $ mpirun -n 128 ./decomposer decomposition128 decomposition 4 8

The above example produces directories ``decomposition512`` (128 * 4) and ``decomposition1024`` (128 * 8).

.. note::
   Parallel decomposition does not re-construct whole problem.
   Decomposition is made on already decomposed parts.
   Hence, multi-level decomposition is **not the same** as simple decomposition,
   even the final number of parts is the same!
   It is strongly recommended to use single decomposition if possible
   in order to avoid load balancing problems.

A decomposed problem can be solved by: ::

  $ mpirun -n 512 ./espreso -i esdata -p decomposition512



OpenFOAM
^^^^^^^^

OpenFOAM support is still under the development.

Problem generator
^^^^^^^^^^^^^^^^^

It is used for a simple generation of a mesh with an arbitrary number of elements and elements types.
This is the key tool for evaluation of scalability of the library at a very large scale.
It is able to generate a multi-billion problem in a few seconds.

Scripts with generator settings can be found in directory ``examples/meshgenerator/``.
Scripts are in the same format as ESPRESO configuration file (``PARAMETER = VALUE``)
and usualy accepts list of nameless command line parameters.
Parameters usualy specify the size of generated problem and element types.
For example, the following command creates and solves a cube composed from
hexahedrons and divided into 8 clusters.
All clusters are divided to 4 subdomains and each subdomain is composed from 16 cubes. ::

  $ mpirun -n 8 ./espreso -p examples/meshgenerator/cube_elasticity_fixed_bottom.txt HEXA8  2 2 2  2 2 1  2 2 4

Check the results
-----------------

The solution can be viewed in `Paraview <http://www.paraview.org/>`__.
The run of ``espreso`` produces output files in the form ``result{MPI_rank}.vtk``.



