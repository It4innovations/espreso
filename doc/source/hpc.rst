

Support for HPC
===============

This section describes scripts for running ESPRESO on clusters.


Solving Ansys Generated Problems on Salomon
-------------------------------------------

To process Ansys generated problems using ESPRESO a section has been added to the following script: ::

  $ ESPRESO_ROOT/machines/salomon.sh

To use this file copy following files to the ESPRESO root directory (the same location where ESPRESO binaries are located): ::

  $ cd ESPRESO_ROOT
  $ cp machines/salomon.sh .
  $ cp machines/ansys_run.conf .

In order to solve a problem following parameters have to be used: ::

  $ ./salomon.sh run_ansys MPI_LIBRARY RUN_TYPE 

The following arguments are used: ::

 run_ansys                 executes Ansys job

 MPI_LIBRARY    intel      use Intel MPI library - recomended 
                sgi        use SGI MPT library	

 RUN_TYPE       mpi        script directly executes mpirun command - usable to interactive sessions 
                pbs        scripts generates PBS jobs and submits a job to PBS


Problem decomposition
---------------------

To process the problem on multiple nodes it has to be decomposed into clusters using the **decomposer** utility: ::
  
  $ ./decomposer input_file destination_directory number_of_clusters1 number_of_clusters2 ..

  Example: 
  $ ./decomposer workbench.dat workbench_decomposed_ 2 4 8

This example will take *workbench.dat* input file and creates three new directories: ::
  
  workbench_decomposed_2
  workbench_decomposed_4
  workbench_decomposed_8

Each of the directory will contain decomposed into particular number of subdomains stored in **ES_DATA** format (internl format of ESPRESO library)

Sample output can looks like this::

   [lriha@login4 espreso]$ ./decomposer /scratch/work/user/lriha/test-small/test_case.dat /scratch/work/user/lriha/test-small/test_case_ 2 4 8
   
   Run ESPRESO on 1 processes.
   Load mesh from Ansys/Workbench format from file /scratch/work/user/lriha/test-small/test_case.dat
   Coordinates loaded - total number of nodes: 729386, average: 729386 (from 729386 to 729386)
   Elements loaded - total number of elements: 393333, average: 393333 (from 393333 to 393333)
   Faces loaded - total number of faces: 0, average: 0 (from 0 to 0)
   Neighbours loaded - number of neighbours for each cluster is average: 0 (from 0 to 0)
   Mesh partitioned into 1 * 1 = 1 parts. There is min 393333, max 393333 elements in subdomain.
   Fix points computed. There is min 0, max 0 fix points in subdomain.
   ********************************************************************************
        Mesh loader         
   ********************************************************************************
   coordinates                                        avg.: 1.137219     min.: 1.137219     max.: 1.137219     sum.: 1.137219     count: 1            % of avg tot: 32.805857   
   elements                                           avg.: 0.971755     min.: 0.971755     max.: 0.971755     sum.: 0.971755     count: 1            % of avg tot: 28.032646   
   faces                                              avg.: 0.000003     min.: 0.000003     max.: 0.000003     sum.: 0.000003     count: 1            % of avg tot: 0.000089    
   boundary conditions                                avg.: 0.005052     min.: 0.005052     max.: 0.005052     sum.: 0.005052     count: 1            % of avg tot: 0.145733    
   cluster boundaries                                 avg.: 0.050290     min.: 0.050290     max.: 0.050290     sum.: 0.050290     count: 1            % of avg tot: 1.450741    
   partition                                          avg.: 1.301362     min.: 1.301362     max.: 1.301362     sum.: 1.301362     count: 1            % of avg tot: 37.540965   
   fix points                                         avg.: 0.000003     min.: 0.000003     max.: 0.000003     sum.: 0.000003     count: 1            % of avg tot: 0.000083    
   corners                                            avg.: 0.000000     min.: 0.000000     max.: 0.000000     sum.: 0.000000     count: 1            % of avg tot: 0.000000    
   --------------------------------------------------------------------------------
   Mesh loader- Total                                 avg.: 3.466512     min.: 3.466512     max.: 3.466512     sum.: 3.466512     count: 1            % of avg tot: 100.000000  
   ********************************************************************************
   Mesh loaded

   Mesh partitiated to 2 parts
   min 358808, max 372704
   Mesh partitiated to 2 parts saved

   Mesh partitiated to 4 parts
   min 176513, max 190326
   Mesh partitiated to 4 parts saved

   Mesh partitiated to 8 parts
   min 82206, max 98083
   Mesh partitiated to 8 parts saved

Basic example to solve the problem
----------------------------------

For this example ESPRESO can be executed as: ::

   $ mpirun -n 4 ./espreso -i esdata -p /scratch/work/user/lriha/test-small/test_case_4 -c ansys5.config -vvv -mmm

Where following parameters are used: ::

  -i esdata                                            specifies input format (in this case ES_DATA)
  -p /scratch/work/user/lriha/test-small/test_case_4   problem to solve with full path to input data
  -c ansys5.config                                     ESPRESO solver configuration file 
  -vvv                                                 full verbose mode 
  -mmm                                                 performs all performance measurements and show results 
 
Please note: This example requires that your environment is setup and you are in a interactive session on a compute node. To avoid this you can use the **salomon.sh** script.

Running the solver using salomon.sh script
------------------------------------------

The **salomon.sh** contains a tool to submit multiple jobs to Salomons PBS queue to measure scalability characteristics. 

This tool has a configuration file **ansys_run.conf**, which is contains following parameters: ::

  ## *****************************************************************************************************************************
  ## Solver directory
  ## *****************************************************************************************************************************

  ## Location where ESPRESO is installed
  ESPRESODIR=~/espreso/



  ## *****************************************************************************************************************************
  ## Input data setup
  ## *****************************************************************************************************************************

  ## Example name without number of clusters
  EXAMPLE=test_case

  ## Example directory location - full path
  EXAMPLE_DIR=/scratch/work/user/lriha/test-small/



  ## *****************************************************************************************************************************
  ## Output data directory setup
  ## *****************************************************************************************************************************

  ## Directory where binaries and setup file are copied and espreso is executed
  WORKDIR=/scratch/work/user/lriha/test-small/results/



  ## *****************************************************************************************************************************
  ## PBS Setup
  ## *****************************************************************************************************************************

  ## PBS queue name
  QUEUE=qmpp

  ## PBS account name
  account=SERVICE



  ## *****************************************************************************************************************************
  ## Expreriment setup
  ## *****************************************************************************************************************************

  ## Number of MPI ranks to be executed - each configuration will be submitted as separate PBS job and stored in separate directory
  ## MPIRANKS=("2" "4" "8" "16" "32" "64" "128" "256")
  MPIRANKS=("2" "4" "8")

  ## Number of MPI processes per compute node
  MPI_PER_NODE=22

  ## Number of threads per MPI process
  THREADS_PER_MPI=1

  ## Solver configuration files - all configuration files will be axecuted in single PBS job
  ## FILES=( "ansys1.config" "ansys2.config" "ansys3.config" "ansys4.config" "ansys5.config")
  FILES=( "ansys.config" )

The script provides two ways to automate the execution of multiple solver runs. This is exposed to user through the two arrays 
in the **ansys_run.conf** file. These are: ::

 - MPIRANKS  - defines the decomposition into clusters. This is primarly used for strong scalability evaluation of the solver.
               Please note: for all configuration an input data has to be prepared with decomposer tool 

 - FILES     - defines different configuration of the solver (different preconditioners, different cluster decomposition, ...) 
               Please note: The different .config files have to be stored in the ESPRESO_ROOT directory 

To provide complete description of the setup, following solver configuration **ansys.config** has been used: ::

  SUBDOMAINS = 32          # each cluster is decomposed into 32 subdomains 
  EPSILON = 1e-2           # stopping criteria in dual  
  FETI_METHOD = 1          # Hybrid Total FETI is used 
  B0_type = 1              # HTFETI corner matrix is composed of face kernels - transformation of basis / averaging is used  
  PRECONDITIONER = 3       # Dirichlet preconditioner 
  ITERATIONS = 10000       # maximum number of iteration 

Now to execute the solver through PBS run: ::

  $ ./salomon.sh run_ansys intel pbs 

The script generates following job submission file for each element of the **MPIRANKS** array in the **ansys_run.conf** and submits it to the PBS queue: :: 

  #!/bin/bash
  export MKL_NUM_THREADS=1
  export OMP_NUM_THREADS=1
  export SOLVER_NUM_THREADS=1
  export PAR_NUM_THREADS=1
  export CILK_NWORKERS=1
  export PARDISOLICMESSAGE=1
  export LD_LIBRARY_PATH=/apps/all/ncurses/5.9-intel-2016.01/lib:/apps/all/zlib/1.2.8-intel-2016.01/lib:/apps/all/tbb/4.4.2.152/tbb/lib/intel64/gcc4.4:/apps/all/imkl/11.3.1.150-iimpi-2016.01-GCC-4.9.3-2.25/mkl/lib/intel64:/apps/all/imkl/11.3.1.150-iimpi-2016.01-GCC-4.9.3-2.25/lib/intel64:/apps/all/impi/5.1.2.150-iccifort-2016.1.150-GCC-4.9.3-2.25/lib64:/apps/all/ifort/2016.1.150-GCC-4.9.3-2.25/compilers_and_libraries_2016.1.150/linux/mpi/intel64:/apps/all/ifort/2016.1.150-GCC-4.9.3-2.25/compilers_and_libraries_2016.1.150/linux/compiler/lib/intel64:/apps/all/ifort/2016.1.150-GCC-4.9.3-2.25/lib/intel64:/apps/all/ifort/2016.1.150-GCC-4.9.3-2.25/lib:/apps/all/icc/2016.1.150-GCC-4.9.3-2.25/compilers_and_libraries_2016.1.150/linux/compiler/lib/intel64:/apps/all/icc/2016.1.150-GCC-4.9.3-2.25/lib/intel64:/apps/all/icc/2016.1.150-GCC-4.9.3-2.25/debugger_2016/libipt/intel64/lib:/apps/all/icc/2016.1.150-GCC-4.9.3-2.25/lib:/apps/all/binutils/2.25-GCCcore-4.9.3/lib:/apps/all/GCCcore/4.9.3/lib/gcc/x86_64-unknown-linux-gnu/4.9.3:/apps/all/GCCcore/4.9.3/lib64:/apps/all/GCCcore/4.9.3/lib:./libs:.

  export LC_CTYPE=

  export MIC_ENV_PREFIX=MIC
  export MIC_OMP_NUM_THREADS=60
  export OFFLOAD_INIT=on_start
  export MIC_USE_2MB_BUFFERS=10k
  export MIC_OMP_NESTED=TRUE
  export MIC_MKL_DYNAMIC=FALSE
  export MIC_MKL_NUM_THREADS=3
  export MIC_OMP_PROC_BIND=spread,close

  module load impi/5.1.2.150-iccifort-2016.1.150-GCC-4.9.3-2.25
  module load imkl/11.3.1.150-iimpi-2016.01-GCC-4.9.3-2.25
  module load tbb/4.4.2.152
  module list

  date | tee -a LOG-test_case_4_22_1_160602_23:23:54.log
  cd /scratch/work/user/lriha/test-small/results//ESP-test_case_4_22_1_160602_23:23:54
  cat $PBS_NODEFILE | tee -a LOG-test_case_4_22_1_160602_23:23:54.node
  mpirun -n 4 ./espreso -i esdata -p /scratch/work/user/lriha/test-small//test_case_4 -c ansys5.config -vvv -mmm | tee -a LOG-test_case_4_22_1_160602_23:23:54.log

For each element of the **MPIRANKS** array in the **ansys_run.conf** a new directory **WORKDIR** (also defined in the **ansys_run.conf**) is created. For this example following directories are created: ::

  $ pwd 

  /scratch/work/user/lriha/test-small/results


  $ ls

  ESP-test_case_2_22_1_160602_23:23:53
  ESP-test_case_4_22_1_160602_23:23:54
  ESP-test_case_8_22_1_160602_23:23:55  

Each directory represent one PBS job and the name can be decoded based on variables from **ansys_run.conf** and date and time as: :: 

  ESP-$EXAMPLE"_"$MPIRANKS[i]"_"$MPI_PER_NODE"_"$THREADS_PER_MPI_$DATE_$TIME

Each directory contains similar content as shows in the following listing: :: 

  $ pwd 

  /scratch/work/user/lriha/test-small/results/ESP-test_case_4_22_1_160602_23:23:54


  $ ls

  # Configuration files of the solver 
  ansys.config
  build.config
  espreso.config

  # Execution scripts 
  salomon.sh
  ansys_run.conf

  # LOG files from PBS 
  es_ans.e794680
  es_ans.o794680

  # ESPRESO log file 
  LOG-test_case_4_22_1_160602_23:23:54.log

  # PBS node file 
  LOG-test_case_4_22_1_160602_23:23:54.node  

  # Job submission script 
  job.qsub

  # ESPRESO Solver binaries and libraries 
  espreso
  libs

  # Solution stored in VTK files 
  result0.vtk
  result1.vtk  
  result2.vtk
  result3.vtk

Results can be opened with Paraview. On Salomon we recomend to use following settings to start the paraview: :: 

  module load ParaView/5.0.0-binary
  module load intel
  #module add CMake/3.0.0-foss-2015g Python/2.7.9-foss-2015g
  export LD_LIBRARY_PATH=/scratch/work/user/sta03/llvm/lib:/scratch/work/user/sta03/openswr-mesa/lib:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/scratch/work/user/sta03/ospray/release:$LD_LIBRARY_PATH
  paraview --mesa-swr-avx2

Then open the **result\*.vtk** files. 

Reducing the size of VTK files
------------------------------

If VTK files are too large to vizualize in Paraview they can be compressed into VTM file using **compVTK.py** script that is part of the ESPRESO. 
This script also reduces the details (by reducing the number surface triangles) of the mesh by almost 90%. 

This script can be executed as: :: 
   
  $ pwd
  /scratch/work/user/lriha/test-small/results/ESP-test_case_4_22_1_160602_23:23:54

  $ ls *.vtk
  result0.vtk  result1.vtk  result2.vtk  result3.vtk
  
  $ mpirun -n 4 pvpython ~/espreso/src/python/compVTK.py 
  size: 4
  rank: 0 to 3
  ['OUT/output3.vtp', 'OUT/output1.vtp', 'OUT/output0.vtp', 'OUT/output2.vtp']

Finally the reduced result can be opened as: :: 

  $ paraview set.vtm 


