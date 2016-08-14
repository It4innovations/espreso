
========================
Support for accelerators
========================


Both the Intel Xeon Phi and GPU acceleration of the computation is supported by the ESPRESO library. 

Intel Xeon Phi
--------------

The ESPRESO library supports an offload of computation to the Intel Xeon Phi coprocessors. The sparse data structures cannot take full advantage of the coprocessor which is equipped with high number (up to 61) SIMD units. To reach the full potential of the architecture a dense representation of the sparse stiffness matrices and their Cholesky decomposition is presented. This is done using local Schur complement matrices obtained from the Pardiso sparse direct solver available in the Intel MKL library.

To enable the support for the accelerators compile and run the library with the following setting. In the `build.config <installation.html#configuration>`_ file set: :: 

  SOLVER = MIC


Change the run time configuration in the ESPRESO `configuration file <run.html#setting-input-parameters>`_ accordingly: ::

  N_MICS = NUMBER_OF_AVAILABLE_MIC_ACCELERATORS
  USE_SCHUR_COMPLEMENT

This ensures that the sparse data structures are replaced by more convenient dense Schur complement matrices and that computation is accelerated using the available coprocessors. Acceleration is supported for both TFETI and HTFETI method. Set the **FETI_METHOD** parameter in the same file accordingly (0 for TFETI, 1 for HTFETI).

To obtain the best performance the following setting of the environmental variables is recommended (see also the `salomon.sh <hpc.html#running-the-solver-using-salomon-sh-script>`_ script available in the installation directory): ::

  export MIC_ENV_PREFIX=MIC
  export MIC_OMP_NUM_THREADS=60
  export MIC_MKL_NUM_THREADS=3
  export MIC_MKL_DYNAMIC=FALSE
  export MIC_OMP_NESTED=TRUE
  export MIC_OMP_PROC_BIND=spread,close
  export MIC_USE_2MB_BUFFERS=10k
  export OFFLOAD_INIT=on_start

If the Intel Xeon Phi acceleration is enabled the local Schur complement are assembled on CPU and offloaded to the coprocessor. If the coprocessor's memory is full the remaining local Schur complements are stored on the host.
  
By default the **load balancing** between a host and coprocessors is enabled. Based on the iteration run time, the program will try to balance workload of the FETI operator application. However, since the arrays storing matrices are duplicated on both host and coprocessors, the approach has higher memory requirements. It is also not suitable for cases when Schur complements do not fit to the coprocessor's memory. To disable load balancing, set ::

  LOAD_BALANCING = 0

in the configuration file.

  
