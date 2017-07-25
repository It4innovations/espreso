
========================
Support for accelerators
========================


Both the Intel Xeon Phi and GPU acceleration of the computation is supported by the ESPRESO library. 

Intel Xeon Phi
--------------

The ESPRESO library supports an offload of computation to the Intel Xeon Phi coprocessors. The sparse data structures cannot take full advantage of the coprocessor which is equipped with high number (up to 61) SIMD units. To reach the full potential of the architecture a dense representation of the sparse stiffness matrices and their Cholesky decomposition is presented. This is done using local Schur complement matrices obtained from the Pardiso sparse direct solver available in the Intel MKL library.

To enable the support for the accelerators compile and run the library with the following setting. In the `build.config <installation.html#configuration>`_ file set: :: 

  SOLVER = MIC

.. warning ::
  ``Salomon`` has incorrectly set MIC_LD_LIBRARY_PATH.
  The correct path is: export MIC_LD_LIBRARY_PATH=/apps/compiler/icc/2016.1.150-GCC-4.9.3/mkl/lib/mic/:/apps/compiler/icc/2016.1.150-GCC-4.9.3/lib/mic/

Change the run time configuration in the ESPRESO `configuration file <run.html#setting-input-parameters>`_ accordingly: ::

  N_MICS = NUMBER_OF_AVAILABLE_MIC_ACCELERATORS
  USE_SCHUR_COMPLEMENT = TRUE

This ensures that the sparse data structures are replaced by more convenient dense Schur complement matrices and that computation is accelerated using the available coprocessors. Acceleration is supported for both TFETI and HTFETI method. Set the **FETI_METHOD** parameter in the same file accordingly (0 for TFETI, 1 for HTFETI).

Alternatively, when setting USE_SCHUR_COMPLEMENT = 0, the MKL version of the Pardiso sparse parallel direct solver will be used on the coprocessor.

To obtain the best performance the following setting of the environmental variables is recommended (see also the `salomon.sh <hpc.html#running-the-solver-using-salomon-sh-script>`_ script available in the installation directory): ::

  export MIC_ENV_PREFIX=MIC
  export MIC_OMP_NUM_THREADS=60
  export MIC_MKL_NUM_THREADS=2
  export MIC_MKL_DYNAMIC=FALSE
  export MIC_OMP_NESTED=TRUE
  export MIC_OMP_PROC_BIND=spread,close
  export MIC_USE_2MB_BUFFERS=10k
  export OFFLOAD_INIT=on_start

If the Intel Xeon Phi acceleration is enabled and the local Schur complement method is used, the local Schur complements are assembled on CPU and offloaded to the coprocessor. If the coprocessor's memory is full the remaining local Schur complements are stored on the host. If the sparse solve is used instead of the local Schur complement method, all subdomain stiffness matrices are factorized on the coprocessor.
  
By default the **load balancing** between a host and coprocessors is enabled. Based on the iteration run time, the program will try to balance workload of the FETI operator application. However, since the arrays storing matrices are duplicated on both host and coprocessors, the approach has higher memory requirements. It is also not suitable for cases when Schur complements do not fit to the coprocessor's memory. To disable load balancing, set ::

  LOAD_BALANCING = FALSE

in the configuration file.

Load balancing is also supported when using the Dirichlet preconditioner. It can be disabled by setting ::

  LOAD_BALANCING_PRECONDITIONER = FALSE

in the configuration file.

Espreso can deal with several configurations of **MPI processes per node**. In the case of one MPI process, this process can utilize all the available coprocessors. In the case when N_MICS = k * ( NUMBER OF MPI PROCESSES PER NODE ), k = 1, 2, ..., each MPI process is assigned k accelerators. In the case when ( NUMBER OF MPI PROCESSES PER NODE ) = l * N_MICS, l = 2, 3, ..., and N_MICS is an even number, each accelerator is shared by ( NUMBER OF MPI PROCESSES PER NODE ) / N_MICS processes. When setting the number of MPI processes per node, keep in mind that within each MPI process an OpenMP thread has to be dedicated per each accelerator to take care of offload.  
