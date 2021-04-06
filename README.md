# ESPRESO 
#### HIGHLY PARALLEL SOLVERS FOR ENGINEERING APPLICATIONS.

ESPRESO is an ExaScale PaRallel FETI SOlver developed at Czech national supercomputing centre IT4Innovations. Main focus of the development team is to create a highly efficient parallel solver which contains several FETI based algorithms including Hybrid Total FETI method suitable for parallel machines with tens or hundreds of thousands of cores. The solver is based on highly efficient communication layer on top of pure MPI.

---
# Instalation
---

###  External Dependencies

In the current version ESPRESO can be compiled and executed on a Linux operating system only. Even the library can be compiled without any other libraries, without [Intel MKL](https://software.intel.com/en-us/intel-mkl) and some of the parallel decomposer, the functionality is strictly reduced. The list of all available wrappers to third party libraries can be queried by:
```sh
$ ./waf -h
```
Libraries should be installed in the system before the ESPRESO installation process. Currently the available wrappers are the followings:
1. Math libraries:
 Currently the ESPRESO core is based on the Intel MKL library. Without this library only MESIO module can be used.
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
 The core the ESPRESO library is based on the in-house FETI based solvers. However, in some circumstances utilization of other solvers can be advantageous (the small number of MPI processes, sequential run, ...). 
 - [Hypre](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods)
 - [Pardiso](https://www.pardiso-project.org/)
5. Other libraries:
 In the case of storing results in the [XDMF](http://www.xdmf.org/index.php/Main_Page) format, HDF5 library should be linked.
 - [HDF5](https://www.hdfgroup.org/solutions/hdf5/)

### Building the library

To compile the ESPRESO library a Python-based framework [Waf](https://waf.io/book/) is used. The compilation process has two phases: **configuration** and **compilation**.  The former configures persistent data and checks all required headers and libraries. The latter builds the library. If something in the environment is changed, it is mandatory to re-run the configuration process.

If all required libraries listed above are available, the following commands build the library:
```sh
$ ./waf configure
$ ./waf
```
The compilation process builds all libraries and executables into the *build* directory. This directory should be added to ``LD_LIBRARY_PATH`` and ``PATH`` environment variables. Then the solver can be started by:
```sh
$ mpirun -n $N espreso -c $ECF
```
where **N** is the number of MPI processes and **ECF** is the ESPRESO configuration file (examples can be found e.g. in the *benchmark* directory).

### Set up the environment

Before running the library, environment variables need to be set. The following variables should be set according to the number of CPU cores per compute node (nCores) and number of MPI processes processed per node (PPN):

 - OMP_NUM_THREADS - should be set to nCores/PPN
 - SOLVER_NUM_THREADS - should be set to nCores/PPN
 - PAR_NUM_THREADS - should be set to nCores/PPN

It is possible to set all environment variables at once usign the following script in the ESPRESO root directory:
```sh
$ . env/threading.default ${nCores/PPN}
```


### Testing the installation
The installation of the library can be validated by the included set of benchmarks which can be executed as follows:
```sh
$ nosetests benchmarks
```
If all tests pass, the library is ready to use.

---
# ESPRESO Solver Interface
---

The library contains a simple C API that allows to utilize FETI based solvers in third party softwares. The API has been successfully tested with an open source multiphysical simulation software [Elmer FEM](https://www.csc.fi/web/elmer)  developed by [CSC - IT Center for Science](https://www.csc.fi/). The interface is available via the feti4i.h header in the include directory.

### Call the library
Usage of the library is described by the simple example. The methods documentation can be found in the library header file.

```cpp
#include "feti4i.h"

int main() {
    // Always initialize MPI before call ESPRESO!!
	MPI_Init(&argc, &argv);

	// Calling of provider should be replaced by appropriate methods in your library !!!
	espreso::APITestESPRESODataProvider provider(&argc, &argv);

	// Solving the problem by FETI4I works in 4 steps:
	//  1. Create the stiffness matrix
	//  2. Configure the ESPRESO library
	//  3. Create an instance of a problem
	//  4. Solve the instance


	// Step 1: Create the stiffness matrix.
	//       : Fill the empty matrix by element data
	FETI4IMatrix K;
	FETI4IInt    matrixType = provider.matrixType();
	FETI4IInt    indexBase = 0;

	FETI4ICreateStiffnessMatrix(&K, matrixType, indexBase); // create matrix
	for (size_t e = 0; e < provider.elements(); e++) {
		provider.addElementMatrix(K, e); // add element data
	}


	// Step 2: Configure the ESPRESO library
	//       : Set all options to default values
	//       : Change required parameters
	FETI4IInt  iopts[FETI4I_INTEGER_OPTIONS_SIZE];
	FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

	FETI4ISetDefaultIntegerOptions(iopts); // set default integer options
	FETI4ISetDefaultRealOptions(ropts); // set default real options

	/* change the default values (see the list of available parameters below)
	iopts[FETI4IIntegerOptions::FETI4I_FETI_METHOD] = 1; // set HYBRID FETI
	iopts[FETI4IIntegerOptions::FETI4I_PRECONDITIONER] = 3; // set Dirichlet preconditioner
	*/


	// Step 3: Create an instance of a problem
	//       : Compute RHS
	//       : Fill L2G, Dirichlet, and neighbours
	FETI4IInstance            instance;
	std::vector<FETI4IReal>   rhs;
	std::vector<FETI4IInt>    dirichlet_indices;
	std::vector<FETI4IReal>   dirichlet_values;
	std::vector<FETI4IInt>    l2g;
	std::vector<FETI4IMPIInt> neighbours;

	provider.computeRHS(rhs);
	provider.fillL2G(l2g);
	provider.fillDirichlet(dirichlet_indices, dirichlet_values);
	provider.fillNeighbours(neighbours);

	FETI4ICreateInstance(
			&instance,
			K, // Stiffness matrix
			rhs.size(), rhs.data(), // RHS
			l2g.data(), // local 2 global indices mapping
			neighbours.size(), neighbours.data(), // neighbours clusters
			dirichlet_indices.size(), dirichlet_indices.data(), dirichlet_values.data(), // Dirichlet boundary condition
			iopts, ropts); // FETI4I options


	// Step 4: Solve the instance
	//       : Process the solution
	std::vector<FETI4IReal> solution(rhs.size());
	FETI4ISolve(instance, solution.size(), solution.data());


	// Finish: Destroy all data
	FETI4IDestroy(K);
	FETI4IDestroy(instance);

	MPI_Finalize();
}
```

### List of available configuration parameters

	FETI4IIntegerOptions::FETI4I_SUBDOMAINS = 4
        The number of subdomains per MPI process

	FETI4IIntegerOptions::FETI4I_MAX_ITERATIONS = 200
	    The maximal number of FETI Solver iterations

	FETI4IIntegerOptions::FETI4I_FETI_METHOD = 0
	    0 = TOTAL_FETI  -- The basis FETI method
	    1 = HYBRID_FETI -- FETI with 2-level decomposition (suitable for huge examples)

	FETI4IIntegerOptions::FETI4I_PRECONDITIONER = 3
	    0 = NONE
	    1 = LUMPED
	    2 = WEIGHT_FUNCTION
	    3 = DIRICHLET

	FETI4IIntegerOptions::FETI4I_CGSOLVER = 0
	    0 = PCG
	    1 = PIPEPCG
	    2 = ORTHOGONALPCG
	    3 = GMRES
	    4 = BICGSTAB

	FETI4IIntegerOptions::FETI4I_N_MICS = 2
	    The number of MIC accelerators per node

	FETI4IIntegerOptions::FETI4I_VERBOSE_LEVEL = 1
	    Verbosity of the library

	FETI4IIntegerOptions::FETI4I_MEASURE_LEVEL = 0
	    Level of time measurement

	FETI4IIntegerOptions::FETI4I_PRINT_MATRICES = 0
	    Turn on/off printing of assembler matrices

    FETI4IRealOptions::FETI4I_PRECISION = 1e-5
        The requested solver precision

### Set default parameters via a configuration file

The default parameters can be changed via a configuration file. The file have to be named `espreso.ecf` and should be placed at the run directory. The structure of the file corresponds with the parameters described above. The next example shows the configuration file with default parameters:

```cpp
FETI4ILIBRARY {
  DOMAINS   4;

  SOLVER {
    #[TOTAL_FETI,HYBRID_FETI]
    METHOD          TOTAL_FETI;
    #[NONE,LUMPED,WEIGHT_FUNCTION,DIRICHLET]
    PRECONDITIONER   DIRICHLET;
    PRECISION            1E-05;
    MAX_ITERATIONS         200;
    #[PCG,PIPEPCG,ORTHOGONALPCG,GMRES,BICGSTAB]
    ITERATIVE_SOLVER       PCG;
    N_MICS                   2;
  }
}
```

---
# MESIO - parallel mesh loader and converter
---

The ESPRESO library contains a highly parallel loader that is able to load unstructured meshes from external sources. The loaded mesh is preprocessed for using with parallel solvers and results are again parallely stored to an external database file (Ensight in default). Besides data preprocessing, this utility can be also used for the parallel conversion among unstructured mesh formats. MESIO is built together with the ESPRESO executable.

Available input unstructured mesh databases:
 - ANSYS CDB
 - OpenFOAM (sequential)
 - VTK Legacy (ASCII)
 - XDMF (mixed data)
 - Ensight (binary format)
  
 Available output unstructured mesh databases:
 - VTK Legacy (ASCII)
 - XDMF (mixed data)
 - Ensight (binary format, default settings)
 - STL surface
 
### Basic configuration via a configuration file

MESIO uses the same configuration files as ESPRESO. However, other then INPUT and OUTPUT parameters are ommited. For both parallel reading the file and parallel decomposition can be set GRANULARITY and REDUCTION_RATIO. These parameters control the number of processes used in the given step (e.g. for reading the input file only by each second MPI process use GRANULARITY=PROCESSES, REDUCTION_RATIO=2; for reading by single MPI process per node use GRANULARITY=NODES, REDUCTION_RATIO=1). 

```cpp
INPUT {
  # [ANSYS_CDB,OPENFOAM,ABAQUS,XDMF,ENSIGHT,VTK_LEGACY]
  FORMAT           ANSYS_CDB; 
  PATH    <path_to_database>;
  
  GRANULARITY      PROCESSES;
  REDUCTION_RATIO          2;
  
  DECOMPOSITION {
    # [PARMETIS,PTSCOTCH]
    PARALLEL_DECOMPOSER   PARMETIS;

    # [NODES,PROCESSES]
    GRANULARITY          PROCESSES;
    REDUCTION_RATIO              1;
  }
}

OUTPUT {
  # [VTK_LEGACY,ENSIGHT,XDMF,STL_SURFACE]
  FORMAT ENSIGHT;
  PATH   results;
}

```
 
# Licence

See the LICENSE file at the root directory.

# Acknowledgement
Development was partially supported by the PRACE 6IP project.
The Partnership for Advanced Computing in Europe (PRACE) is an international non-profit association with its seat in Brussels.
The PRACE Research Infrastructure provides a persistent world-class high performance computing service for scientists and researchers from academia and industry in Europe.
The computer systems and their operations accessible through PRACE are provided by 5 PRACE members (BSC representing Spain, CINECA representing Italy, ETH Zurich/CSCS representing Switzerland, GCS representing Germany and GENCI representing France).
The Implementation Phase of PRACE receives funding from the EU’s Horizon 2020 Research and Innovation Programme (2014-2020) under grant agreement 823767. For more information, see www.prace-ri.eu.
