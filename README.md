The main repository of the highly-parallel framework for engineering applications developed at IT4Innovations National Supercomputing Center (http://numbox.it4i.cz/). Main focus of the development team is to create a highly efficient parallel tools suitable for parallel machines with tens or hundreds of thousands of cores. The added value of the framework are our
 - mesh pre/post-processing tool [mesio](#mesio) with ability to load and store various sequential unstructured mesh formats in parallel, and
 - massively parallel linear solver [espreso](#espreso) based on FETI Domain Decomposition Methods, which is capable of solving problems with over a hundred billions of unknowns using thousands compute nodes.

# MESIO
Mesio is a highly-parallel loader and convertor of external unstructured meshes databases. It is composed of several commonly used algorithms that together provide a robust solution for the pre/post-processing phase of solving large-scale engineering problems ([pdf](https://doi.org/10.1109/IPDPS.2019.00084)). With mesio a user is able to use the same sequential database file with an arbitrary number of MPI processes. It allows engineers to use their favourite tools for creation of numerical models without any penalty for a parallel run.

Mesio is composed from (i) the mesh builder that is able to reconstruct an unstructured mesh from randomly scattered data across parallel processes in a distributed memory of an HPC machine without gathering data into a single node and (ii) lightweight parallel parsers of given mesh databases. Hence, the usage is not limited to any particular format -- only a simple parallel parser has to be implemented for a new format. Currently, the following parsers are implemented:
 - ANSYS CDB
 - Ensight
 - VTK Legacy
 - XDMF
 - Netgen
 - OpenFOAM (partially)
 - Abaqus (partially)

An output database stored by mesio is also in a sequetial form for simple by a favourite visualization tool. The following format are available:
 - VTK Legacy
 - XDMF
 - Ensight
 - STL surface

# ESPRESO
Espreso is a set of several highly-scalable solvers based on the methods of domain decomposition designed to take full advantage of today's most powerful petascale supercomputers ([pdf](https://dx.doi.org/10.1177/1094342018798452)). It contains in-house developed FETI based algorithms including the Hybrid Total FETI method suitable for parallel machines with tens or hundreds of thousands of cores.

The solver also contains a general [API](#espreso-api) for usage in a third party software. The API has been successfully used to connect ESPRESO with the Elmer ([pdf](https://dx.doi.org/10.1007/978-3-319-97136-0_10)). Even though ESPRESO is C++ library, the API uses plain C only. Hence, it is easy to use it with various other languages such as Fortran.

Both mesio and espreso can be configured by [ecf](#configuration-via-a-configuration-file) files. In addition, espreso can be also configured by a simple [GUI](#espreso-gui).

---
---
---

## Instalation
---

####  External Dependencies

In the current version the modules can be compiled and executed on Linux operating systems only. Some functionality requires third party libraries that should be installed in the system before the installation process. Currently available wrappers are the followings:
1. Math libraries:
 Currently the espreso core is based on the Intel MKL library. Without this library only mesio module can be used.
 * [Intel MKL](https://software.intel.com/en-us/intel-mkl) (at least version 2018.4)
2. Parallel graph partitioners:
 For loading external databases the library should be linked with an external parallel graph partitioner.
 - [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
 - [PT-Scotch](https://www.labri.fr/perso/pelegrin/scotch/)
3. Sequential graph partitioners:
 In other to use second level decomposition (e.q. Hybrid FETI) at external graph decomposer is needed. Except KaHIP, they are usually part of parallel graph decomposers.
 - [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
 - [Scotch](https://www.labri.fr/perso/pelegrin/scotch/)
 - [KaHIP](http://algo2.iti.kit.edu/kahip/)
4. Third party solvers:
 The espreso is based on in-house FETI based solvers. However, in some circumstances utilization of other solvers can be advantageous (the small number of MPI processes, sequential run, etc.).
 - [Hypre](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods)
 - [Pardiso](https://www.pardiso-project.org/)
5. Other libraries:
 In the case of loading or storing results in the [XDMF](http://www.xdmf.org/index.php/Main_Page) format, HDF5 library should be linked. Qt is needed for building the GUI module.
 - [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
 - [Qt5](https://www.qt.io/)

#### Building the library

For compilation Python-based framework [Waf](https://waf.io/book/) is used. The compilation process has two phases: **configuration** and **compilation**.  The former configures persistent data and checks available headers and libraries. The latter builds the library. It is mandatory to re-run the configuration process after any environment change. The following commands build all modules if require libraries are available:
```sh
$ ./waf configure
$ ./waf
```
The compilation process builds all libraries and executables into the *build* directory. This directory should be added to ``LD_LIBRARY_PATH`` and ``PATH`` environment variables. Then it is possible to run [mesio](#mesio) or [espreso](#espreso) by the following command:
```sh
$ mpirun -n $N mesio -c $ECF
$ mpirun -n $N espreso -c $ECF
```
where **N** is the number of MPI processes and **ECF** is the ESPRESO configuration file (examples can be found e.g. in the *benchmarks* directory).

#### Set up the environment

Before running the library, the following variables should be set according to the number of CPU cores per compute node (nCores) and number of MPI processes processed per node (PPN):

 - OMP_NUM_THREADS - should be set to nCores/PPN
 - SOLVER_NUM_THREADS - should be set to nCores/PPN
 - PAR_NUM_THREADS - should be set to nCores/PPN

It is possible to set all environment variables at once usign the following script in the ESPRESO root directory:
```sh
$ . env/threading.default ${nCores/PPN}
```


#### Testing the installation
The installation of the library can be validated by the included set of benchmarks. Usually it is sufficient to run *dummy* benchmarks only in order to avoid exhausting testing. The testing by [nosetests](https://nose.readthedocs.io/en/latest/) can be started by the following command:
```sh
$ nosetests benchmarks/dummy
```

---
---
---

## Usage of the tools
---

#### Configuration via a configuration file

Run of libraries is controlled by a configuration file. By this file a user is able to manage every single parameter of the libraries. The following file is an example of the configuration file for loading a mesh from an ANSYS CDB database and set the solver for heat transfer (see *benchmarks* directory for more examples).
```cpp
INPUT {
  FORMAT           ANSYS_CDB;
  PATH    <path_to_database>;
}

OUTPUT {
  FORMAT ENSIGHT;
  PATH   results;
}

PHYSICS   HEAT_TRANSFER_3D;
HEAT_TRANSFER_3D {
  LOAD_STEPS        1;

  MATERIALS {
    MATERIAL_1 {
      DENS 1;  CP 1;

      THERMAL_CONDUCTIVITY {
        MODEL   DIAGONAL;
        KXX  1; KYY 10; KZZ 10;
      }
    }
  }

  MATERIAL_SET {
    ALL_ELEMENTS   MATERIAL_1;
  }

  INITIAL_TEMPERATURE {
    ALL_ELEMENTS   200;
  }

  LOAD_STEPS_SETTINGS {
    1 {
      DURATION_TIME     1;
      TYPE   STEADY_STATE;
      MODE         LINEAR;
      SOLVER         FETI;

      FETI {
        METHOD          TOTAL_FETI;
        PRECONDITIONER   DIRICHLET;
        PRECISION            1E-08;
        ITERATIVE_SOLVER       PCG;
      }

      TEMPERATURE {
        TOP 100; BOTTOM 300;
      }
    }
  }
}
```

#### ESPRESO API

The library contains a simple C API that allows to utilize FETI based solvers in third party softwares. The interface is available via the `feti4i.h` header in the include directory and the `libfeti4i` library. An example of calling the library can be found in `src/api/apitester.cpp`. The tester uses the espreso physical module to assemble matrices and call the FETI solver through API.

Usage of the library is scatched by the following code. The methods documentation can be found in the library header file.

```cpp
#include "feti4i.h"

int main(int argc, char **argv)
{
	/* -----------------------*/
	// initialize the library //
	/* -----------------------*/
	MPI_Init(&argc, &argv);
	FETI4IInit(MPI_COMM_WORLD, 0); // VERBOSE_LEVEL = 0

	/* ---------------------------------------------------------------------------------------------------------*/
	// set options (see 'src/config/ecf/linearsolver/feti.h' for all available options of the given parameters) //
	/* ---------------------------------------------------------------------------------------------------------*/
	FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	FETI4ISetDefaultIntegerOptions(iopts);
	FETI4ISetDefaultRealOptions(ropts);

	iopts[FETI4I_SUBDOMAINS] = 8;
	iopts[FETI4I_MAX_ITERATIONS] = 1000;
	iopts[FETI4I_FETI_METHOD] = 0; // 0: TOTAL FETI; 1: HYBRID FETI
	iopts[FETI4I_PRECONDITIONER] = 3; // 0: NONE; 1: LUMPED; 2: WEIGHT; 3: DIRICHLET
	iopts[FETI4I_CGSOLVER] = 0; // 0: PCG; 1: Pipelined PCG; 2: orthogonal PCG; 3: GMRES; 4: BICGSTAB
	ropts[FETI4I_PRECISION] = 1e-8;

    /* ---------------------------------------------------------------------*/
	// create stiffness matrix and provide local-to-global mapping of nodes //
	/* ---------------------------------------------------------------------*/
	FETI4IMatrix matrix;
	FETI4ICreateStiffnessMatrix(&matrix, nnodes, l2g, indexing, matrix_type, dofs_per_node, FETI4I_ORDER_GROUPED_DOFS);
	for (int e = 0; e < elements; ++e) {
	    FETI4IAddElement(matrix, FETI4I_ETYPE_VOLUME, 8, nodes, stiffness);
	}

	/* ------------------------------------------------------------------------*/
	// create instance of the system (it also call FETI solver pre-processing) //
	/* ------------------------------------------------------------------------*/
	FETI4IInstance            instance;
	FETI4ICreateInstance(&instance, matrix, neighbors_size, neighbors, dirichlet_size, dirichlet_indices dirichlet_values, iopts, ropts);

	/* ---------------------------------*/
	// compute RHS and solve the system //
	/* ---------------------------------*/
	FETI4ISolve(instance, rhs, solution);

	FETI4IDestroy(matrix);
	FETI4IDestroy(instance);
	FETI4IFinalize();

	MPI_Finalize();
	return 0;
}
```

#### ESPRESO GUI

It is possible to use GUI, instead of setting a configuration file manually. GUI is automatically generated from internal configuration structures that assures always up to date list of available parameters. Besides showing all parameters with their possible values, GUI also showes the input geometry, and elements and boundary regions. It is also possible to use GUI for generation of a configuration file.

In order to built GUI, one has to configure espreso with flag *\-\-with-gui*. In that case Qt5 has to be installed.

```sh
$ ./waf configure --with-gui
$ ./waf
```

# Licence

See the LICENSE file at the root directory.

# Acknowledgement

This work was supported by
 - The Ministry of Education, Youth and Sports from the Large Infrastructures for Research, Experimental Development and Innovations project "e-Infrastructure CZ -- LM2018140",
 - The Ministry of Education, Youth and Sports from the National Programme of Sustainability (NPS II) project "IT4Innovations excellence in science -- LQ1602",
 - The IT4Innovations infrastructure which is supported from the Large Infrastructures for Research, Experimental Development and Innovations project "IT4Innovations National Supercomputing Center -- LM2015070".

Development was also partially supported by the PRACE 6IP project. The Partnership for Advanced Computing in Europe (PRACE) is an international non-profit association with its seat in Brussels. The PRACE Research Infrastructure provides a persistent world-class high performance computing service for scientists and researchers from academia and industry in Europe. The computer systems and their operations accessible through PRACE are provided by 5 PRACE members (BSC representing Spain, CINECA representing Italy, ETH Zurich/CSCS representing Switzerland, GCS representing Germany and GENCI representing France). The Implementation Phase of PRACE receives funding from the EU’s Horizon 2020 Research and Innovation Programme (2014-2020) under grant agreement 823767. For more information, see www.prace-ri.eu.

