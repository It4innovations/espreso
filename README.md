
The main repository of the highly-parallel framework for engineering applications developed at IT4Innovations National Supercomputing Center (http://numbox.it4i.cz/).
Main focus of the development team is to create a highly efficient parallel tools suitable for parallel machines with tens or hundreds of thousands of cores. The added value of the framework are our:
 - Massively parallel linear solver [espreso](#espreso) based on FETI Domain Decomposition Methods, which is capable of solving problems with over a hundred billions of unknowns using thousands compute nodes.
 - Mesh pre/post-processing tool [mesio](#mesio) with ability to load and store various sequential unstructured mesh formats in parallel.

# ESPRESO
Espreso is a set of several highly-scalable solvers based on the methods of domain decomposition designed to take full advantage of today's most powerful petascale supercomputers ([pdf](https://dx.doi.org/10.1177/1094342018798452)). It contains in-house developed FETI based algorithms including the Hybrid Total FETI method suitable for parallel machines with tens or hundreds of thousands of cores. The solver also supports GPU acceleration by [CUDA](https://developer.nvidia.com/cuda-toolkit), [ROCm](https://www.amd.com/en/products/software/rocm.html), and [oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.edajzv).



# MESIO
Mesio is a highly-parallel loader and converter of external unstructured meshes databases. It provides an efficient solution for the pre/post-processing phase of solving large-scale engineering problems ([pdf](https://doi.org/10.1109/IPDPS.2019.00084)). With mesio a user is able to use the same sequential database file with an arbitrary number of MPI processes. It allows engineers to use their favorite tools for creation of numerical models without any penalty for a parallel run.

---
---
---

## Dependencies

Before the installation the library, please install the following third party libraries:

 1. Math libraries are fundamental for the espreso library. Install at least one of the following option:

  - [Intel MKL](https://software.intel.com/en-us/intel-mkl)
    - at least version 2019 is required
  - [BLAS](https://www.netlib.org/blas/), [LAPACK](https://www.netlib.org/lapack/), [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html)
    - set 'BLAS_LIBRARIES' to installed BLAS libraries (e.g., BLAS_LIBRARIES=openblas)
    - set 'LAPACK_LIBRARIES' to installed LAPACK libraries (e.g., LAPACK_LIBRARIES=openblas)
    - at least version 7.6.0 of SuiteSparse is required

 2. Graph partitioning tool are fundamental for the mesio library. Install at least one of the following option:

  - [METIS](https://github.com/KarypisLab/METIS), [ParMETIS](https://github.com/KarypisLab/ParMETIS)
  - [Scotch](https://www.labri.fr/perso/pelegrin/scotch/), [PT-Scotch](https://www.labri.fr/perso/pelegrin/scotch/)

 3. Distributed sparse solvers can be used instead of in-house FETI solver for. Note that they are suitable for small examples only as the scalability of these solvers are limited.

  - [MKL Parallel Direct Sparse Solver](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/onemkl-pardiso-parallel-direct-sparse-solver-iface.html)

 4. Espreso supports GPU acceleration with all mayor GPU vendors. Once found in the system, particular version is automatically installed.

  - [CUDA](https://developer.nvidia.com/cuda-toolkit)
    - at least version 11.7.0
  - [ROCm](https://www.amd.com/en/products/software/rocm.html)
    - at least version 5.4.3
  - [oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html#gs.edajzv)
    - at least version 2024.1

 5. Multi-physical coupling is allowed with the following library:

  - [Precice](https://precice.org/)

## Installation

For compilation Python-based framework [Waf](https://waf.io/book/) is used. The compilation process has two phases: **configuration** and **compilation**.  The former configures persistent data and checks available of third party libraries. The latter builds the library. It is mandatory to re-run the configuration process after any environment change. Run the following commands to build the library:
```sh
$ ./waf configure
$ ./waf
```
The compilation process builds all libraries and executable tools into the *build* directory. This directory should be added to ``LD_LIBRARY_PATH`` and ``PATH`` environment variables. Then it is possible to run [mesio](#mesio) or [espreso](#espreso) by the following command:
```sh
$ mpirun -n N mesio -c ECF
$ mpirun -n N espreso -c ECF
```
where **N** is the number of MPI processes and **ECF** is the ESPRESO configuration file (examples can be found e.g. in the *tests* directory).

## Testing the installation
The installation of the library can be validated by the included set of tests. The testing by [nose2](https://docs.nose2.io/en/latest/) can be started by the following command:
```sh
$ nose2 -v -s tests/
```

## Configuration via a configuration file

Run of libraries is controlled by a configuration file. By this file a user is able to manage every single parameter of the libraries. The following file is an example of the configuration file for loading a mesh from an ANSYS CDB database and set the solver for heat transfer (see *tests* directory for more examples).
```cpp
INPUT {
  FORMAT           ANSYS_CDB;
  PATH    <path_to_database>;
}

OUTPUT {
  FORMAT ENSIGHT;
  PATH   results;
}

PHYSICS   HEAT_TRANSFER;
HEAT_TRANSFER {
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
## MESIO API

Our parallel loader can be utilized by third party software by provided C API. The interface is available via the `mesio.h` header in the include directory and the `libmesioapi` library. An example of calling the library can be found in `src/api/api.mesio.cpp`. The code below shows how the loader should be called. Once the method `MESIOLoad` is finished an input database is loaded. Then, one can use provided functions to return mesh data stored in internal structures. The code shows requesting of nodes and elements only. For the full API description see the provided example.

```cpp
#include "mesio.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

	MESIO mesio;
	MESIOInit(MPI_COMM_WORLD, 2);
	MESIOLoad(&mesio, MESIO_ANSYS, "path_to_file", MESIO_PARMETIS, 4);

	{ // GET NODES
		MESIOInt nhalo, offset, size, totalSize;
		MESIOInt *ids, *position;
		MESIOReal *coordinates;

		MESIONodes(mesio, &nhalo, &offset, &size, &totalSize, &ids, &position, &coordinates);
	}

	{ // GET ELEMENTS
		MESIOInt offset, size, totalSize;
		MESIOInt *type, *enodesDist, *enodesData;

		MESIOElements(mesio, &offset, &size, &totalSize, &type, &enodesDist, &enodesData);
	}
	MESIOFinalize();
	MPI_Finalize();
	return 0;
}
```

## ESPRESO API

Our FETI based solvers can be utilized by third party software by provided C API. The interface is available via the `feti4i.h` header in the include directory and the `libfeti4i` library. An example of calling the library can be found in `src/api/api.feti4i.cpp`. The tester uses the espreso physical module to assemble matrices and call the FETI solver through API.

Usage of the library is scratched by the following code. The methods documentation can be found in the library header file.

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

# License

See the LICENSE file at the root directory.

# Acknowledgment

This work was supported by
 - The Ministry of Education, Youth and Sports from the Large Infrastructures for Research, Experimental Development and Innovations project "e-Infrastructure CZ -- LM2018140",
 - The Ministry of Education, Youth and Sports from the National Programme of Sustainability (NPS II) project "IT4Innovations excellence in science -- LQ1602",
 - The IT4Innovations infrastructure which is supported from the Large Infrastructures for Research, Experimental Development and Innovations project "IT4Innovations National Supercomputing Center -- LM2015070".

Development was also partially supported by the PRACE 6IP project. The Partnership for Advanced Computing in Europe (PRACE) is an international non-profit association with its seat in Brussels. The PRACE Research Infrastructure provides a persistent world-class high performance computing service for scientists and researchers from academia and industry in Europe. The computer systems and their operations accessible through PRACE are provided by 5 PRACE members (BSC representing Spain, CINECA representing Italy, ETH Zurich/CSCS representing Switzerland, GCS representing Germany and GENCI representing France). The Implementation Phase of PRACE receives funding from the EU’s Horizon 2020 Research and Innovation Programme (2014-2020) under grant agreement 823767. For more information, see www.prace-ri.eu.
