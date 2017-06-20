

ESPRESO API
===========

To interface the ESPRESO with the third party software, a general API has been implemented.
The API has been successfully used to connect the ESPRESO with the `Elmer <https://csc.fi/web/elmer/elmer>`_.

Even though the ESPRESO is C++ library, the API uses plain C only.
Hence, it is easy to use it with various other languages such as Fortran.

To setup and solve problem using the API is done in the following steps: 

 - create a stiffness matrix,
 - create an instance of a problem,
 - configure the solver,
 - run the solver.

A simple working example can be found at ``libespreso/apiexample.cpp``.
Fo simplicity the example is written in C++:

.. code-block:: cpp
   :linenos:
   :emphasize-lines: 22, 24, 28, 31 - 35, 38 - 41, 44 - 57, 63, 67 - 68

   // Always initialize the MPI before calling the ESPRESO
   MPI_Init(&argc, &argv);

   // The following data are considered to be computed by a library
   std::vector<std::vector<FETI4IInt> >  K_indices;
   std::vector<std::vector<FETI4IReal> > K_values;
   std::vector<FETI4IReal>               rhs;
   std::vector<FETI4IInt>                dirichlet_indices;
   std::vector<FETI4IReal>               dirichlet_values;
   std::vector<FETI4IInt>                l2g;
   std::vector<FETI4IMPIInt>             neighbours;

   // The input data can be loaded from the ESPRESO example directory
   loadStructures(
      "../examples/api/cube",
      K_indices, K_values, rhs,
      dirichlet_indices, dirichlet_values,
      l2g, neighbours);


   // Step 1: Create the stiffness matrix
   FETI4IMatrix K;
   FETI4IInt indexBase = 0;
   FETI4ICreateStiffnessMatrix(&K, indexBase);

   // Compose the matrix from elements matrices
   for (size_t i = 0; i < K_indices.size(); i++) {
      FETI4IAddElement(K, K_indices[i].size(), K_indices[i].data(), K_values[i].data());
   }

   // Step 2: Configure the ESPRESO solver
   FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
   FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

   FETI4ISetDefaultIntegerOptions(iopts);
   FETI4ISetDefaultRealOptions(ropts);

   iopts[FETI4I_SUBDOMAINS] = 8;
   iopts[FETI4I_PRECONDITIONER] = 3;
   iopts[FETI4I_VERBOSE_LEVEL] = 3;
   iopts[FETI4I_MEASURE_LEVEL] = 3;

   // Step 3: Create an instance of a problem
   FETI4IInstance instance;
   FETI4ICreateInstance(
         &instance,
         K,
         rhs.size(),
         rhs.data(),
         l2g.data(),
         neighbours.size(),
         neighbours.data(),
         dirichlet_indices.size(),
         dirichlet_indices.data(),
         dirichlet_values.data(),
         iopts,
         ropts);


   // Step 4: Solve the problem
   // Allocate the memory for the solution
   std::vector<FETI4IReal> solution(rhs.size());

   // Solve the system
   FETI4ISolve(instance, solution.size(), solution.data());

   // Process the solution
   // ... 

   FETI4IDestroy(K);
   FETI4IDestroy(instance);

   MPI_Finalize();

The highlighted lines are the ESPRESO API functions.
At first, a stiffness matrix has to be created (lines 22, 23).
The method ``FETI4ICreateStiffnessMatrix`` (line 24)
accepts the data holder ``FETI4IMatrix`` and ``indexBase``.
When the stiffness matrix data holder is created,
the element matrices can be added by the ``FETI4IAddElement`` method (line 28).

The ESPRESO settings are defined by an array of integers and an array of floating point values (lines 31, 32).
These arrays have to be passed to the ``FETI4ICreateInstance`` method. 
User should setup the ESPRESO to the default values by calling the ``FETI4ISetDefaultIntegerOptions`` and ``FETI4ISetDefaultRealOptions`` functions.
Then it can change any of the parameters to the required value.


The stiffness matrix and other required data structures are passed to the ESPRESO solver
by calling the method ``FETI4ICreateInstance`` (lines 33 - 43). The method returns a data holder for the created instance.

.. note::

   The instance is prepared according to the passed parameters. Hence, a later change of
   parameters or options has no effect on the instance.
   Instance can be changed by the API ``Update`` method.

The problem can be solved by the ``FETI4ISolve`` method (line 49).
The solution is then saved to the vector that needs to be preallocated (line 46).

The data holders should be destroyed by calling the ``FETI4IDestroy`` method.


.. note::

   **neighbours** is an array of neighbouring MPI ranks. Hence, the first rank is 0.

   **dirichlet_indices** are in the local numbering.

A description of API methods
----------------------------

.. doxygenfunction:: FETI4ICreateStiffnessMatrix

.. doxygenfunction:: FETI4IAddElement

.. doxygenfunction:: FETI4ICreateInstance

.. doxygenfunction:: FETI4ISolve

.. doxygenfunction:: FETI4IDestroy

