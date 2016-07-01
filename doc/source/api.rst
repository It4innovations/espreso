

ESPRESO API
===========

For usage with other software, a general API was developed.
The API was successfully tested with cooperation with `Elmer <https://csc.fi/web/elmer/elmer>`_.
Even the ESPRESO is C++ library, the API use only plain C.
Hence, it is simple to use it with various other languages (e.g. Fortran).

Usage of the API contains three steps:
create stiffness matrix,
create an instance of a problem,
and run the solver.
A simple working example can be found at ``libespreso/apiexample.cpp``.
The example is written in C++ for its simplicity:

.. code-block:: cpp
   :linenos:
   :emphasize-lines: 22, 24, 28, 31 - 35, 38 - 41, 44 - 57, 63, 67 - 68

   // Always initialize MPI before call ESPRESO!
   MPI_Init(&argc, &argv);

   // The following data are considered to be computed by a library
   std::vector<std::vector<FETI4IInt> >  K_indices;
   std::vector<std::vector<FETI4IReal> > K_values;
   std::vector<FETI4IReal>               rhs;
   std::vector<FETI4IInt>                dirichlet_indices;
   std::vector<FETI4IReal>               dirichlet_values;
   std::vector<FETI4IInt>                l2g;
   std::vector<FETI4IMPIInt>             neighbours;

   // We load data stored in ESPRESO example directory
   loadStructures(
      "../examples/api/cube",
      K_indices, K_values, rhs,
      dirichlet_indices, dirichlet_values,
      l2g, neighbours);


   // At first create stiffness matrix
   FETI4IMatrix K;
   FETI4IInt indexBase = 0;
   FETI4ICreateStiffnessMatrix(&K, indexBase);

   // Compose the matrix from elements matrices
   for (size_t i = 0; i < K_indices.size(); i++) {
      FETI4IAddElement(K, K_indices[i].size(), K_indices[i].data(), K_values[i].data());
   }

   FETI4IInt iopts[FETI4I_INTEGER_OPTIONS_SIZE];
   FETI4IReal ropts[FETI4I_REAL_OPTIONS_SIZE];

   FETI4ISetDefaultIntegerOptions(iopts);
   FETI4ISetDefaultRealOptions(ropts);

   // Configure ESPRESO solver
   iopts[FETI4I_SUBDOMAINS] = 8;
   iopts[FETI4I_PRECONDITIONER] = 3;
   iopts[FETI4I_VERBOSE_LEVEL] = 3;
   iopts[FETI4I_MEASURE_LEVEL] = 3;

   // Create instance of a problem
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

   // Prepare memory for save solution
   std::vector<FETI4IReal> solution(rhs.size());

   // Solve the system
   FETI4ISolve(instance, solution.size(), solution.data());

   // Process solution

   FETI4IDestroy(K);
   FETI4IDestroy(instance);

   MPI_Finalize();

The highlighted lines are ESPRESO API functions.
At first, a stiffness matrix has to be created (lines 22, 23).
The method ``FETI4ICreateStiffnessMatrix`` (line 24)
accepts a data holder ``FETI4IMatrix`` and ``indexBase``.
When the stiffness matrix data holder is created,
element matrices can be added by method ``FETI4IAddElement`` (line 28).

The settings of ESPRESO is controlled by arrays of integer and real values (lines 31, 32).
These arrays have to be passed to ``FETI4ICreateInstance``. You can set all parameters
to appropriate values (see `complete list <parameters.html>`__ for the description of parameters)
or you can use default values set by ``FETI4ISetDefaultIntegerOptions`` (``FETI4ISetDefaultRealOptions``).

A filled stiffness matrix and other needed data can passed to the ESPRESO solver
by method ``FETI4ICreateInstance`` (lines 33 - 43).
Again, the method returns a data holder for a created instance.

.. note::

   The instance is prepared according to passed parameters. Hence, a later change of
   parameters or options has no effect on the instance.
   Instance can be changed by API ``Update`` methods.

The instance can be solved by method ``FETI4ISolve`` (line 49).
The solution is saved to prepared vector (line 46).
Data holders should be destroyed by ``FETI4IDestroy``.

.. note::

   **neighbours** is an array of neighbours MPI ranks. Hence, the first rank is 0.

   **dirichlet_indices** are in local numbering.

A description of API methods
----------------------------

.. doxygenfunction:: FETI4ICreateStiffnessMatrix

.. doxygenfunction:: FETI4IAddElement

.. doxygenfunction:: FETI4ICreateInstance

.. doxygenfunction:: FETI4ISolve

.. doxygenfunction:: FETI4IDestroy

