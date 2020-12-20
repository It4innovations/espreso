#include "superclusteracc.h"

#include "assembler/instance.h"

using namespace espreso;

SuperClusterAcc::~SuperClusterAcc() {
    if (this->deleteMatrices) {
        for (esint i = 0; i < configuration.n_mics; i++) {
            if (matricesPerAcc[i]) {
                //delete [] matricesPerAcc[i];
            }
        }
    }
}


void SuperClusterAcc::SetAcceleratorAffinity() {

    //        	Logical to Physical Processor Mapping
    //        	? Hardware:
    //        	? Physical Cores are 0..60
    //        	? Logical Cores are 0..243
    //        	? Mapping is not what you are used to!
    //        	? Logical Core 0 maps to Physical core 60, thread context 0
    //        	? Logical Core 1 maps to Physical core 0, thread context 0
    //        	? Logical Core 2 maps to Physical core 0, thread context 1
    //        	? Logical Core 3 maps to Physical core 0, thread context 2
    //        	? Logical Core 4 maps to Physical core 0, thread context 3
    //        	? Logical Core 5 maps to Physical core 1, thread context 0
    //        	? ...
    //        	? Logical Core 240 maps to Physical core 59, thread context 3
    //        	? Logical Core 241 maps to Physical core 60, thread context 1
    //        	? Logical Core 242 maps to Physical core 60, thread context 2
    //        	? Logical Core 243 maps to Physical core 60, thread context 3
    //        	? OpenMP threads start binding to logical core 1, not logical core 0
    //        	? For compact mapping 240 OpenMP threads are mapped to the first 60 cores
    //        	? No contention for the core containing logical core 0 ? the core that the O/S uses most
    //        	? But for scatter and balanced mappings, contention for logical core 0 begins at 61 threads
    //        	? Not much performance impact unless O/S is very busy
    //        	? Best to avoid core 60 for offload jobs & MPI jobs with compute/communication overlap
    //        	? KMP_PLACE_THREADS limits the range of cores & thread contexts
    //        	? E.g., KMP_PLACE_THREADS=60c,2c with KMP_AFFINITY=compact and OMP_NUM_THREADS=120 places 2
    //        	threads on each of the first 60 cores


    // detect how many MPI processes is running per node
    int _MPIglobalRank;
    int _MPIglobalSize;
    int _MPInodeRank;
    int _MPInodeSize;
    std::string _nodeName;
    MPI_Comm _currentNode;
    MPI_Comm _storingProcs;
    MPI_Comm_rank(info::mpi::comm, &_MPIglobalRank);
    MPI_Comm_size(info::mpi::comm, &_MPIglobalSize);
    int color;
    int length;
    std::vector<char> name(MPI_MAX_PROCESSOR_NAME);
    MPI_Get_processor_name(name.data(), &length);
    _nodeName = std::string(name.begin(), name.begin() + length);
    std::vector<int> rCounts(_MPIglobalSize);
    std::vector<int> rDispl(_MPIglobalSize);
    std::vector<int> colors(_MPIglobalSize);
    std::vector<char> names;
    MPI_Gather(&length, sizeof(int), MPI_BYTE, rCounts.data(), sizeof(int), MPI_BYTE, 0, info::mpi::comm);
    for (size_t i = 1; i < rCounts.size(); i++) {
        rDispl[i] += rDispl[i - 1] + rCounts[i - 1];
    }
    names.resize(rDispl.back() + rCounts.back());
    MPI_Gatherv(name.data(), length * sizeof(char), MPI_BYTE, names.data(), rCounts.data(), rDispl.data(), MPI_BYTE, 0, info::mpi::comm);
    std::map<std::string, size_t> nodes;
    for (size_t i = 0; i < _MPIglobalSize; i++) {
        std::string str(names.begin() + rDispl[i], names.begin() + rDispl[i] + rCounts[i]);
        auto it = nodes.find(str);
        if (it == nodes.end()) {
            size_t s = nodes.size();
            nodes[str] = s;
        }
        colors[i] = nodes[str];
    }
    MPI_Scatter(colors.data(), sizeof(int), MPI_BYTE, &color, sizeof(int), MPI_BYTE, 0, info::mpi::comm);
    MPI_Comm_split(info::mpi::comm, color, _MPIglobalRank, &_currentNode);
    MPI_Comm_rank(_currentNode, &_MPInodeRank);
    MPI_Comm_size(_currentNode, &_MPInodeSize);
    MPI_Comm_split(info::mpi::comm, _MPInodeRank, _MPIglobalRank, &_storingProcs);

    this->MPI_per_node = _MPInodeSize;
    // END - detect how many MPI processes is running per node
    ESINFO(PROGRESS2) << "MPI ranks per node: " << _MPInodeSize;

    int nMICs = configuration.n_mics;

    if ( _MPInodeSize > nMICs && nMICs % 2 == 0 && _MPInodeSize % nMICs == 0 ) {
        // the case when there is more MPInodeSize = 2*k*nMICs
        this->MPI_per_acc = _MPInodeSize / nMICs;
        this->acc_per_MPI = 1;
        this->myTargets.reserve( 1 );
        this->myTargets.push_back( _MPInodeRank / MPI_per_acc );
        this->acc_rank = _MPInodeRank % this->MPI_per_acc;
    } else if ( nMICs % _MPInodeSize == 0 ) {
        // the case when 2*k*MPInodeSize = nMICS
        this->MPI_per_acc = 1;
        this->acc_per_MPI = nMICs / _MPInodeSize;
        this->myTargets.reserve( this->acc_per_MPI );
        for ( esint i = 0; i < acc_per_MPI ; ++i ) {
            this->myTargets.push_back( _MPInodeRank * acc_per_MPI + i );    
        }
        this->acc_rank = 0;
    } else if ( nMICs == 1 && _MPInodeSize > 1 ) {
        this->MPI_per_acc = _MPInodeSize;
        this->acc_per_MPI = 1;
        this->myTargets.reserve( 1 );
        this->myTargets.push_back( 0 );
        this->acc_rank = _MPInodeRank;
    } else {
        ESINFO(PROGRESS2) << "Incorrect number of MPI processes per accelerator!" << _MPInodeSize;  
    }

#pragma omp parallel num_threads( acc_per_MPI )
    {
        int used_core_num = 0;
        int first_core = 0;
        int last_core = 0;
        int target = myTargets.at(omp_get_thread_num());
        int MPIperAcc = this->MPI_per_acc;
        int rank = this->acc_rank;

#pragma offload target(mic:target) 
        {
            cpu_set_t my_set;        // Define cpu_set bit mask. 
            CPU_ZERO(&my_set);       // Initialize it all to 0
            int cores = sysconf(_SC_NPROCESSORS_ONLN); // for Xeon Phi 7120 - results is 244
            cores = ( cores / 4 ) - 1; // keep core for system and remove effect of hyperthreading
            int cores_per_rank = cores / MPIperAcc;

            for (int i = 0; i < cores_per_rank ; i++) {
                if (i == 0) {
                    //first_core = 1*(cores_per_rank)*rank + 1*i;
                    first_core = cores_per_rank*rank + 1;
                }
                last_core = 1*(cores_per_rank)*rank + 1*i;
                //int core = 1 + 4*(cores_per_rank)*rank + 4*i;

                int core =1+  4*cores_per_rank*rank + 4*i;
                for (int j = 0 ; j < 4; j++ ) {
                    CPU_SET(core , &my_set);     /* set the bit that represents core 7. */
                    core++;
                }
                used_core_num++;
            }   

            sched_setaffinity(0, sizeof(cpu_set_t), &my_set); /* Set affinity of tihs process to */
            omp_set_num_threads(4*cores_per_rank);
            /* the defined mask, i.e. only 7. */
        }
    }

}


void SuperClusterAcc::Create_SC_perDomain(bool USE_FLOAT) {
    // Ratio of work done on MIC

    ESINFO(PROGRESS3) << "Creating Local Schur complements";
    double MICr = 1.0;
    if ( configuration.load_balancing ) {
        MICr = 0.9;
    }
    // First, get the available memory on coprocessors (in bytes)
    double usableRAM = 0.8;
    long *micMem = new long[ this->acc_per_MPI ];

    int target = 0;
    for (esint i = 0; i < this->acc_per_MPI; ++i) {
        long currentMem = 0;
        target = this->myTargets.at(i);
#pragma offload target(mic:target)
        {
            long pages = sysconf(_SC_AVPHYS_PAGES);
            long page_size = sysconf(_SC_PAGE_SIZE);
            currentMem = pages * page_size;
        }
        micMem[i] = currentMem / this->MPI_per_acc;
    }


#pragma omp parallel for
    for (esint i = 0; i < number_of_subdomains_per_supercluster; i++ ) {
        domains[i]->B1_comp_dom.MatTranspose(domains[i]->B1t_comp_dom);
    }

    // compute sizes of data to be offloaded to MIC
    esint *matrixPerPack = new esint[this->acc_per_MPI];

    for (esint i = 0 ; i < this->acc_per_MPI; ++i) {
        matrixPerPack[i] = domains.size() / this->acc_per_MPI;
    }
    for (esint i = 0 ; i < domains.size() % this->acc_per_MPI; ++i) {
        matrixPerPack[i]++;
    }
    esint offset = 0;
    bool symmetric = SYMMETRIC_SYSTEM;
    this->B1KplusPacks.resize( this->acc_per_MPI, DenseMatrixPack(configuration) );
    this->accDomains.resize( this->acc_per_MPI );

    esint scalarSize = (USE_FLOAT) ? sizeof(float) : sizeof(double);

    for ( int i = 0; i < this->acc_per_MPI; i++ )
    {
        long dataSize = 0;
        long currDataSize = 0;
        bool MICFull = false;
        esint numCPUDomains = 0;

        for ( int j = offset; j < offset + matrixPerPack[i]; j++ ) {
            if (!symmetric) {
                currDataSize = domains[j]->B1t_comp_dom.cols * domains[j]->B1t_comp_dom.cols;
            } else {
                // isPacked => is symmetric
                currDataSize = ( ( 1.0 + ( double ) domains[j]->B1t_comp_dom.cols ) *
                        ( ( double ) domains[j]->B1t_comp_dom.cols ) / 2.0 );
            }
            long dataInBytes = (currDataSize + dataSize) * scalarSize;
            if (MICFull || (dataInBytes > usableRAM * micMem[i])) {
                // when no more memory is available at MIC leave domain on CPU
                MICFull = true;
                hostDomains.push_back(j);
                numCPUDomains++;
            } else {
                accDomains[i].push_back(j);
                dataSize += currDataSize;
            }
        }

        // it is necessary to subtract numCPUDomains AFTER setting offset!
        offset += matrixPerPack[i];
        matrixPerPack[i] -= numCPUDomains;

        this->B1KplusPacks[i].Resize( matrixPerPack[i], dataSize, USE_FLOAT );
        this->B1KplusPacks[i].setMICratio( MICr );

        if ( configuration.load_balancing ) {
            this->B1KplusPacks[i].enableLoadBalancing();
        } else {
            this->B1KplusPacks[i].disableLoadBalancing();
        }

        for ( esint j = 0; j < matrixPerPack[i]; ++j ) {
            this->B1KplusPacks[ i ].PreparePack( j, domains[accDomains[i].at(j)]->B1t_comp_dom.cols,
                    domains[accDomains[i].at(j)]->B1t_comp_dom.cols, symmetric );
        }
    }

    delete [] micMem;

    // iterate over available MICs and assemble their domains on CPU
    for (esint i = 0 ; i < this->acc_per_MPI ; ++i) {
        target = this->myTargets.at(i);
        this->B1KplusPacks[i].AllocateVectors( );
        this->B1KplusPacks[i].SetDevice( target );
#pragma omp parallel for
        for (esint j = 0 ; j < accDomains[i].size() ; ++j) {
            esint domN = accDomains[i].at(j);
            SparseMatrix tmpB;
            domains[domN]->B1_comp_dom.MatTranspose(tmpB);
            SparseSolverCPU tmpsps;
            tmpsps.Create_SC_w_Mat( domains[domN]->K, tmpB, domains[domN]->B1Kplus, false, symmetric);
            if (!USE_FLOAT) {

                double * matrixPointer = this->B1KplusPacks[i].getMatrixPointer(j);
                memcpy(matrixPointer, &(domains[domN]->B1Kplus.dense_values[0]),
                        this->B1KplusPacks[i].getDataLength(j) * sizeof(double) );
                SEQ_VECTOR<double>().swap(  domains[domN]->B1Kplus.dense_values);
            } else {

                float * matrixPointer = this->B1KplusPacks[i].getMatrixPointer_fl(j) ;
                domains[domN]->B1Kplus.ConvertDenseToDenseFloat( 1 );
                //for (esint k = 0; k << this->B1KplusPacks[i].getDataLength(j); ++k) {
                //    matrixPointer[k] = domains[domN].B1Kplus.dense_values_fl[k];                 
                //    }
                memcpy(matrixPointer, &(domains[domN]->B1Kplus.dense_values_fl[0]),
                        this->B1KplusPacks[i].getDataLength(j) * sizeof(float) );
                SEQ_VECTOR<float>().swap(  domains[domN]->B1Kplus.dense_values_fl);

            }
            ESINFO(PROGRESS3) << Info::plain() << ".";
        }
    }

#pragma omp parallel num_threads(this->acc_per_MPI)
    {
        this->B1KplusPacks[omp_get_thread_num()].CopyToMIC();
    }

#pragma omp parallel for
    for (esint d = 0; d < hostDomains.size(); ++d) {

        ESINFO(PROGRESS3) << Info::plain() << "*";
        esint domN = hostDomains.at(d);
        SparseMatrix TmpB;
        domains[domN]->B1_comp_dom.MatTranspose(TmpB);
        SparseSolverCPU tmpsps;
        tmpsps.Create_SC_w_Mat( domains[domN]->K, TmpB, domains[domN]->B1Kplus, false, symmetric);
        if (USE_FLOAT) {
            domains[domN]->B1Kplus.ConvertDenseToDenseFloat(1);
            domains[domN]->B1Kplus.USE_FLOAT = true;
        }
    }


    delete [] matrixPerPack;

}


void SuperClusterAcc::SetupKsolvers ( ) {
    // this part is for setting CPU pardiso, temporarily until everything is
    // solved on MIC
#pragma omp parallel for
    for (esint d = 0; d < domains.size(); d++) {

        // Import of Regularized matrix K into Kplus (Sparse Solver)
        switch (configuration.Ksolver) {
            case FETIConfiguration::KSOLVER::DIRECT_DP:
                domains[d]->Kplus.ImportMatrix_wo_Copy (domains[d]->K);
                break;
            case FETIConfiguration::KSOLVER::ITERATIVE:
                domains[d]->Kplus.ImportMatrix_wo_Copy (domains[d]->K);
                break;
            case FETIConfiguration::KSOLVER::DIRECT_SP:
                domains[d]->Kplus.ImportMatrix_fl(domains[d]->K);
                break;
            case FETIConfiguration::KSOLVER::DIRECT_MP:
                domains[d]->Kplus.ImportMatrix_fl(domains[d]->K);
                break;
                //		case 4:
                //			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
                //			break;
            default:
                ESINFO(ERROR) << "Invalid KSOLVER value.";
                exit(EXIT_FAILURE);
        }

        if (configuration.keep_factors) {
            std::stringstream ss;
            ss << "init -> rank: " << info::mpi::rank << ", subdomain: " << d;
            domains[d]->Kplus.keep_factors = true;
            if (configuration.Ksolver != FETIConfiguration::KSOLVER::ITERATIVE) {
                domains[d]->Kplus.Factorization (ss.str());
            }
        } else {
            domains[d]->Kplus.keep_factors = false;
            domains[d]->Kplus.MPIrank = info::mpi::rank;
        }

        domains[d]->domain_prim_size = domains[d]->Kplus.cols;
        //TODO: Hot Fix - needs to be done better
        if ( !SYMMETRIC_SYSTEM ) {
            // 11 = Real and unsymmetric matrix
            domains[d]->Kplus.mtype = 11;
        } else {
            // 2 = Real and symmetric positive definite
            domains[d]->Kplus.mtype = 2;
        }
        //TODO: else stokes = -2 = Real and symmetric indefinite

        if ( d == 0 && info::mpi::rank == 0) {
            //    domains[d]->Kplus.msglvl = Info::report(LIBRARIES) ? 1 : 0;
        }
        domains[d]->Kplus.msglvl = 0;
    }
    if (!USE_KINV) {
        double MICr = 0.8;

        esint *matrixPerPack = new esint[ this->acc_per_MPI ];

        for ( esint i = 0; i < this->acc_per_MPI; ++i ) {
            matrixPerPack[i] = domains.size() / this->acc_per_MPI;    
        }
        for ( esint i = 0 ; i < domains.size() % this->acc_per_MPI; ++i ) {
            matrixPerPack[i]++;
        }

        esint offset = 0; 
        bool use_float = (configuration.Ksolver == FETIConfiguration::KSOLVER::DIRECT_SP);

        this->SparseKPack.resize( this->acc_per_MPI, SparseMatrixPack(configuration, use_float) );
        this->accDomains.resize( this->acc_per_MPI );
        this->matricesPerAcc.reserve( this->acc_per_MPI );

        for ( int i = 0; i < this->acc_per_MPI; ++i ) {
            this->matricesPerAcc[ i ] = new SparseMatrix*[ matrixPerPack[ i ] ];
            for ( int j = offset; j < offset + matrixPerPack[ i ]; ++j ) {
                accDomains[ i ].push_back( j );       
                this->matricesPerAcc[ i ][ j - offset ] = &(domains[j]->K);
            }

            offset += matrixPerPack[ i ];

            this->SparseKPack[ i ].AddMatrices( this->matricesPerAcc[ i ], 
                    matrixPerPack[ i ], this->myTargets.at( i ) );
            this->SparseKPack[ i ].setMICratio( MICr );

            if ( configuration.load_balancing ) {
                this->SparseKPack[ i ].enableLoadBalancing( );
            } else {
                this->SparseKPack[ i ].disableLoadBalancing( );
            }
        }

#pragma omp parallel num_threads( acc_per_MPI )
        {
            this->SparseKPack[ omp_get_thread_num() ].AllocateVectors();
            this->SparseKPack[ omp_get_thread_num() ].CopyToMIC();
            this->SparseKPack[ omp_get_thread_num() ].FactorizeMIC();
        }
        delete [] matrixPerPack;

    }
}


void SuperClusterAcc::SetupPreconditioner() {

    if ( configuration.preconditioner != FETIConfiguration::PRECONDITIONER::DIRICHLET ) {
        for ( esint c = 0; c < clusters.size(); ++c ) {
            clusters[ c ].SetupPreconditioner();
        }
    } else {
        CreateDirichletPrec( instance ); 
    }
}


void SuperClusterAcc::CreateDirichletPrec( DataHolder *instance ) {

    bool USE_FLOAT = ( configuration.schur_precision == FETIConfiguration::FLOAT_PRECISION::SINGLE ||
            configuration.Ksolver == FETIConfiguration::KSOLVER::DIRECT_SP );

    // Ratio of work done on MIC
    double MICr = 1.0;
    if ( configuration.load_balancing_preconditioner ) {
        MICr = 0.9;
    }

    // First, get the available memory on coprocessors (in bytes)
    double usableRAM = 0.8;
    long *micMem = new long[this->acc_per_MPI];

    int target = 0;
    for (esint i = 0; i < this->acc_per_MPI; ++i) {
        long currentMem = 0;
        target = this->myTargets.at(i);
#pragma offload target(mic:target)
        {
            long pages = sysconf(_SC_AVPHYS_PAGES);
            long page_size = sysconf(_SC_PAGE_SIZE);
            currentMem = pages * page_size;
        }
        micMem[i] = currentMem / this->MPI_per_acc;
    }

    // compute sizes of data to be offloaded to MIC
    esint *matrixPerPack = new esint[this->acc_per_MPI];

    for (esint i = 0 ; i < this->acc_per_MPI; ++i) {
        matrixPerPack[i] = domains.size() / this->acc_per_MPI;
    }
    for (esint i = 0 ; i < domains.size() % this->acc_per_MPI; ++i) {
        matrixPerPack[i]++;
    }
    esint offset = 0;
    bool symmetric = SYMMETRIC_SYSTEM;
    this->DirichletPacks.resize( this->acc_per_MPI, DenseMatrixPack(configuration) );
    this->accPreconditioners.resize( this->acc_per_MPI );

    esint scalarSize = (USE_FLOAT) ? sizeof(float) : sizeof(double);

    for ( int i = 0; i < this->acc_per_MPI; i++ )
    {
        long dataSize = 0;
        long currDataSize = 0;
        bool MICFull = false;
        esint numCPUDomains = 0;

        for ( int j = offset; j < offset + matrixPerPack[i]; j++ ) {
            esint size = domains[j]->B1t_Dir_perm_vec.size();
            if (!symmetric) {
                currDataSize = size * size;
            } else {
                // isPacked => is symmetric
                currDataSize = ( ( 1.0 + ( double ) size ) *
                        ( ( double ) size ) / 2.0 );
            }
            long dataInBytes = (currDataSize + dataSize) * scalarSize;
            if (MICFull || (dataInBytes > usableRAM * micMem[i])) {
                // when no more memory is available at MIC leave domain on CPU
                MICFull = true;
                hostPreconditioners.push_back(j);
                numCPUDomains++;
            } else {
                accPreconditioners[i].push_back(j);
                dataSize += currDataSize;
            }
        }

        // it is necessary to subtract numCPUDomains AFTER setting offset!
        offset += matrixPerPack[i];
        matrixPerPack[i] -= numCPUDomains;

        this->DirichletPacks[i].Resize( matrixPerPack[i], dataSize, USE_FLOAT );
        this->DirichletPacks[i].setMICratio( MICr );

        if ( configuration.load_balancing_preconditioner ) {
            this->DirichletPacks[i].enableLoadBalancing();
        } else {
            this->DirichletPacks[i].disableLoadBalancing();
        }

        for ( esint j = 0; j < matrixPerPack[i]; ++j ) {
            this->DirichletPacks[ i ].PreparePack( j, domains[accPreconditioners[i].at(j)]->B1t_Dir_perm_vec.size(),
                    domains[accPreconditioners[i].at(j)]->B1t_Dir_perm_vec.size(), symmetric );
        }
    }

    delete [] micMem;
    delete [] matrixPerPack;

    for (esint mic = 0 ; mic < this->acc_per_MPI ; ++mic ) {
        target = this->myTargets.at(mic);
        this->DirichletPacks[mic].AllocateVectors( );
        this->DirichletPacks[mic].SetDevice( target );        

#pragma omp parallel for
        for (esint j = 0; j < accPreconditioners[mic].size(); ++j ) {
            esint d = accPreconditioners[mic].at(j);
            SEQ_VECTOR <esint> perm_vec = domains[d]->B1t_Dir_perm_vec;
            SEQ_VECTOR <esint> perm_vec_full ( instance->K[d].rows );
            SEQ_VECTOR <esint> perm_vec_diff ( instance->K[d].rows );

            SEQ_VECTOR <esint> I_row_indices_p (instance->K[d].nnz);
            SEQ_VECTOR <esint> J_col_indices_p (instance->K[d].nnz);

            for (esint i = 0; i < perm_vec.size(); i++) {
                perm_vec[i] = perm_vec[i] - 1;
            }

            for (esint i = 0; i < perm_vec_full.size(); i++) {
                perm_vec_full[i] = i;
            }

            auto it = std::set_difference( perm_vec_full.begin(), perm_vec_full.end(), perm_vec.begin(), perm_vec.end(), perm_vec_diff.begin() );
            perm_vec_diff.resize(it - perm_vec_diff.begin());

            perm_vec_full = perm_vec_diff;
            perm_vec_full.insert(perm_vec_full.end(), perm_vec.begin(), perm_vec.end());

            SparseMatrix K_modif = instance->K[d];
            SparseMatrix RegMatCRS = instance->RegMat[d];
            RegMatCRS.ConvertToCSRwithSort(0);
            K_modif.MatAddInPlace(RegMatCRS,'N',-1);
            // K_modif.RemoveLower();

            SEQ_VECTOR <SEQ_VECTOR<esint >> vec_I1_i2(K_modif.rows, SEQ_VECTOR<esint >(2, 1));
            esint offset = K_modif.CSR_I_row_indices[0] ? 1 : 0;

            for (esint i = 0; i < K_modif.rows;i++){
                vec_I1_i2[i][0] = perm_vec_full[i];
                vec_I1_i2[i][1] = i; // position to create reverse permutation
            }

            std::sort(vec_I1_i2.begin(), vec_I1_i2.end(), [](const SEQ_VECTOR <esint >& a, const SEQ_VECTOR<esint>& b) { return a[0] < b[0]; });

            // permutations made on matrix in COO format
            K_modif.ConvertToCOO(0);
            esint I_index,J_index;
            bool unsymmetric=!SYMMETRIC_SYSTEM;
            for (esint i = 0;i<K_modif.nnz;i++){
                I_index = vec_I1_i2[K_modif.I_row_indices[i]-offset][1]+offset;
                J_index = vec_I1_i2[K_modif.J_col_indices[i]-offset][1]+offset;
                if (unsymmetric || I_index<=J_index){
                    I_row_indices_p[i]=I_index;
                    J_col_indices_p[i]=J_index;
                }
                else{
                    I_row_indices_p[i]=J_index;
                    J_col_indices_p[i]=I_index;
                }
            }
            for (esint i = 0; i<K_modif.nnz;i++){
                K_modif.I_row_indices[i] = I_row_indices_p[i];
                K_modif.J_col_indices[i] = J_col_indices_p[i];
            }
            K_modif.ConvertToCSRwithSort(1);
            {
                if (info::ecf->env.print_matrices) {
                    std::ofstream osS(Logging::prepareFile(d, "K_modif"));
                    osS << K_modif;
                    osS.close();
                }
            }


            // ------------------------------------------------------------------------------------------------------------------
            bool diagonalized_K_rr = configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET;
            //        PRECONDITIONER==NONE              - 0
            //        PRECONDITIONER==LUMPED            - 1
            //        PRECONDITIONER==WEIGHT_FUNCTION   - 2
            //        PRECONDITIONER==DIRICHLET         - 3
            //        PRECONDITIONER==SUPER_DIRICHLET   - 4
            //        
            //        When next line is uncomment, var. PRECONDITIONER==DIRICHLET and PRECONDITIONER==SUPER_DIRICHLET provide identical preconditioner.
            //        bool diagonalized_K_rr = false
            // ------------------------------------------------------------------------------------------------------------------

            esint sc_size = perm_vec.size();

            if (sc_size == instance->K[d].rows) {
                domains[d]->Prec = instance->K[d];
                domains[d]->Prec.ConvertCSRToDense(1);
                // if physics.K[d] does not contain inner DOF
            } else {

                if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET || 
                        configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET ) {
                    SparseSolverCPU createSchur;
                    //          createSchur.msglvl=1;
                    esint sc_size = perm_vec.size();
                    createSchur.ImportMatrix_wo_Copy(K_modif);
                    createSchur.Create_SC(domains[d]->Prec, sc_size,false);
                    domains[d]->Prec.ConvertCSRToDense(1);
                }
                else
                {
                    SparseMatrix K_rr;
                    SparseMatrix K_rs;
                    SparseMatrix K_sr;
                    SparseMatrix KsrInvKrrKrs; 

                    esint i_start = 0;
                    esint nonsing_size = K_modif.rows - sc_size - i_start;
                    esint j_start = nonsing_size;

                    K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, nonsing_size,j_start,sc_size);

                    if (SYMMETRIC_SYSTEM){
                        K_rs.MatTranspose(K_sr);
                    }
                    else
                    {
                        K_sr.getSubBlockmatrix_rs(K_modif,K_sr,j_start,sc_size,i_start, nonsing_size);
                    }

                    domains[d]->Prec.getSubDiagBlockmatrix(K_modif,domains[d]->Prec,nonsing_size,sc_size);
                    SEQ_VECTOR <double> diagonals;
                    SparseSolverCPU K_rr_solver;

                    // K_rs is replaced by:
                    // a) K_rs = 1/diag(K_rr) * K_rs          (simplified Dirichlet precond.)
                    // b) K_rs =    inv(K_rr) * K_rs          (classical Dirichlet precond. assembled by own - not via PardisoSC routine)
                    if (diagonalized_K_rr){
                        diagonals = K_modif.getDiagonal();
                        // diagonals is obtained directly from K_modif (not from K_rr to avoid assembling) thanks to its structure
                        //      K_modif = [K_rr, K_rs]
                        //                [K_sr, K_ss]
                        // 
                        for (esint i = 0; i < K_rs.rows; i++) {
                            for (esint j = K_rs.CSR_I_row_indices[i]; j < K_rs.CSR_I_row_indices[i + 1]; j++) {
                                K_rs.CSR_V_values[j - offset] /= diagonals[i];
                            }
                        }
                    }
                    else
                    {
                        K_rr.getSubDiagBlockmatrix(K_modif,K_rr,i_start, nonsing_size);
                        K_rr_solver.ImportMatrix_wo_Copy(K_rr);
                        //            K_rr_solver.msglvl = 1;
                        K_rr_solver.SolveMat_Dense(K_rs);
                    }

                    KsrInvKrrKrs.MatMat(K_sr,'N',K_rs);
                    domains[d]->Prec.MatAddInPlace(KsrInvKrrKrs,'N',-1);
                    //          if (!diagonalized_K_rr){
                    //				    domains[d].Prec.ConvertCSRToDense(1);
                    //          }
                }

            }

            if (info::ecf->env.print_matrices) {
                std::ofstream osS(Logging::prepareFile(d, "S"));
                SparseMatrix SC =  domains[d]->Prec;
                if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET || 
                        configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET ){
                    SC.ConvertDenseToCSR(1);
                }
                osS << SC;
                osS.close();
            }

            // insert matrix to MatrixPack and delete it
            if (!USE_FLOAT) {
                double * matrixPointer = this->DirichletPacks[mic].getMatrixPointer( j );
                memcpy( matrixPointer, &( domains[ d ]->Prec.dense_values[ 0 ] ), 
                        this->DirichletPacks[ mic ].getDataLength( j ) * sizeof( double ) );
                SEQ_VECTOR<double>().swap( domains[d]->Prec.dense_values );
            } else {
                float * matrixPointer = this->DirichletPacks[mic].getMatrixPointer_fl( j );
                domains[d]->Prec.ConvertDenseToDenseFloat( 1 );
                memcpy( matrixPointer, &(domains[d]->Prec.dense_values_fl[0]), this->DirichletPacks[ mic ].getDataLength(j) * sizeof(float) );
            }
            ESINFO(PROGRESS3) << Info::plain() << ".";
        }
    }

#pragma omp parallel num_threads(this->acc_per_MPI)
    {
        this->DirichletPacks[omp_get_thread_num()].CopyToMIC();
    }


    // finish the blocks held on CPU
#pragma omp parallel for
    for (esint j = 0; j < hostPreconditioners.size(); ++j ) {
        esint d = hostPreconditioners.at(j);
        SEQ_VECTOR <esint> perm_vec = domains[d]->B1t_Dir_perm_vec;
        SEQ_VECTOR <esint> perm_vec_full ( instance->K[d].rows );
        SEQ_VECTOR <esint> perm_vec_diff ( instance->K[d].rows );

        SEQ_VECTOR <esint> I_row_indices_p (instance->K[d].nnz);
        SEQ_VECTOR <esint> J_col_indices_p (instance->K[d].nnz);

        for (esint i = 0; i < perm_vec.size(); i++) {
            perm_vec[i] = perm_vec[i] - 1;
        }

        for (esint i = 0; i < perm_vec_full.size(); i++) {
            perm_vec_full[i] = i;
        }

        auto it = std::set_difference( perm_vec_full.begin(), perm_vec_full.end(), perm_vec.begin(), perm_vec.end(), perm_vec_diff.begin() );
        perm_vec_diff.resize(it - perm_vec_diff.begin());

        perm_vec_full = perm_vec_diff;
        perm_vec_full.insert(perm_vec_full.end(), perm_vec.begin(), perm_vec.end());

        SparseMatrix K_modif = instance->K[d];
        SparseMatrix RegMatCRS = instance->RegMat[d];
        RegMatCRS.ConvertToCSRwithSort(0);
        K_modif.MatAddInPlace(RegMatCRS,'N',-1);
        // K_modif.RemoveLower();

        SEQ_VECTOR <SEQ_VECTOR<esint >> vec_I1_i2(K_modif.rows, SEQ_VECTOR<esint >(2, 1));
        esint offset = K_modif.CSR_I_row_indices[0] ? 1 : 0;

        for (esint i = 0; i < K_modif.rows;i++){
            vec_I1_i2[i][0] = perm_vec_full[i];
            vec_I1_i2[i][1] = i; // position to create reverse permutation
        }

        std::sort(vec_I1_i2.begin(), vec_I1_i2.end(), [](const SEQ_VECTOR <esint >& a, const SEQ_VECTOR<esint>& b) { return a[0] < b[0]; });

        // permutations made on matrix in COO format
        K_modif.ConvertToCOO(0);
        esint I_index,J_index;
        bool unsymmetric=!SYMMETRIC_SYSTEM;
        for (esint i = 0;i<K_modif.nnz;i++){
            I_index = vec_I1_i2[K_modif.I_row_indices[i]-offset][1]+offset;
            J_index = vec_I1_i2[K_modif.J_col_indices[i]-offset][1]+offset;
            if (unsymmetric || I_index<=J_index){
                I_row_indices_p[i]=I_index;
                J_col_indices_p[i]=J_index;
            }
            else{
                I_row_indices_p[i]=J_index;
                J_col_indices_p[i]=I_index;
            }
        }
        for (esint i = 0; i<K_modif.nnz;i++){
            K_modif.I_row_indices[i] = I_row_indices_p[i];
            K_modif.J_col_indices[i] = J_col_indices_p[i];
        }
        K_modif.ConvertToCSRwithSort(1);
        {
            if (info::ecf->env.print_matrices) {
                std::ofstream osS(Logging::prepareFile(d, "K_modif"));
                osS << K_modif;
                osS.close();
            }
        }


        // ------------------------------------------------------------------------------------------------------------------
        bool diagonalized_K_rr = configuration.preconditioner == FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET;
        //        PRECONDITIONER==NONE              - 0
        //        PRECONDITIONER==LUMPED            - 1
        //        PRECONDITIONER==WEIGHT_FUNCTION   - 2
        //        PRECONDITIONER==DIRICHLET         - 3
        //        PRECONDITIONER==SUPER_DIRICHLET   - 4
        //        
        //        When next line is uncomment, var. PRECONDITIONER==DIRICHLET and PRECONDITIONER==SUPER_DIRICHLET provide identical preconditioner.
        //        bool diagonalized_K_rr = false
        // ------------------------------------------------------------------------------------------------------------------

        esint sc_size = perm_vec.size();

        if (sc_size == instance->K[d].rows) {
            domains[d]->Prec = instance->K[d];
            domains[d]->Prec.ConvertCSRToDense(1);
            // if physics.K[d] does not contain inner DOF
        } else {

            if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET) {
                SparseSolverCPU createSchur;
                //          createSchur.msglvl=1;
                esint sc_size = perm_vec.size();
                createSchur.ImportMatrix_wo_Copy(K_modif);
                createSchur.Create_SC(domains[d]->Prec, sc_size,false);
                domains[d]->Prec.ConvertCSRToDense(1);
            }
            else
            {
                SparseMatrix K_rr;
                SparseMatrix K_rs;
                SparseMatrix K_sr;
                SparseMatrix KsrInvKrrKrs; 

                esint i_start = 0;
                esint nonsing_size = K_modif.rows - sc_size - i_start;
                esint j_start = nonsing_size;

                K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, nonsing_size,j_start,sc_size);

                if (SYMMETRIC_SYSTEM){
                    K_rs.MatTranspose(K_sr);
                }
                else
                {
                    K_sr.getSubBlockmatrix_rs(K_modif,K_sr,j_start,sc_size,i_start, nonsing_size);
                }

                domains[d]->Prec.getSubDiagBlockmatrix(K_modif,domains[d]->Prec,nonsing_size,sc_size);
                SEQ_VECTOR <double> diagonals;
                SparseSolverCPU K_rr_solver;

                // K_rs is replaced by:
                // a) K_rs = 1/diag(K_rr) * K_rs          (simplified Dirichlet precond.)
                // b) K_rs =    inv(K_rr) * K_rs          (classical Dirichlet precond. assembled by own - not via PardisoSC routine)
                if (diagonalized_K_rr){
                    diagonals = K_modif.getDiagonal();
                    // diagonals is obtained directly from K_modif (not from K_rr to avoid assembling) thanks to its structure
                    //      K_modif = [K_rr, K_rs]
                    //                [K_sr, K_ss]
                    // 
                    for (esint i = 0; i < K_rs.rows; i++) {
                        for (esint j = K_rs.CSR_I_row_indices[i]; j < K_rs.CSR_I_row_indices[i + 1]; j++) {
                            K_rs.CSR_V_values[j - offset] /= diagonals[i];
                        }
                    }
                }
                else
                {
                    K_rr.getSubDiagBlockmatrix(K_modif,K_rr,i_start, nonsing_size);
                    K_rr_solver.ImportMatrix_wo_Copy(K_rr);
                    //            K_rr_solver.msglvl = 1;
                    K_rr_solver.SolveMat_Dense(K_rs);
                }

                KsrInvKrrKrs.MatMat(K_sr,'N',K_rs);
                domains[d]->Prec.MatAddInPlace(KsrInvKrrKrs,'N',-1);
                //          if (!diagonalized_K_rr){
                //				    domains[d]->Prec.ConvertCSRToDense(1);
                //          }
            }

        }

        if (info::ecf->env.print_matrices) {
            std::ofstream osS(Logging::prepareFile(d, "S"));
            SparseMatrix SC =  domains[d]->Prec;
            if (configuration.preconditioner == FETIConfiguration::PRECONDITIONER::DIRICHLET){
                SC.ConvertDenseToCSR(1);
            }
            osS << SC;
            osS.close();
        }

        if (USE_FLOAT) {
            domains[d]->Prec.ConvertDenseToDenseFloat( 1 );
            domains[d]->Prec.USE_FLOAT = true;
        }

        ESINFO(PROGRESS3) << Info::plain() << ".";
    }

    ESINFO(PROGRESS3);   
}


void SuperClusterAcc::multKplusGlobal_l_Acc(SEQ_VECTOR<SEQ_VECTOR<double> *> & x_in, 
        double & CPUtime, 
        double * MICtime ) {


    SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster;
    x_prim_cluster.resize(number_of_subdomains_per_supercluster);

		for (int c = 0; c < numClusters; c++) {

			for (int d = 0; d < clusters[c].domains.size(); d++)
				x_prim_cluster[d].swap( *x_in[clusters[c].domains[d].domain_global_index] );

            clusters[c].multKplusGlobal_l_prepare_Acc( x_prim_cluster );

			for (int d = 0; d < clusters[c].domains.size(); d++)
				(*x_in[clusters[c].domains[d].domain_global_index]).swap( x_prim_cluster[d] );

		}


    esint maxDevNumber = this->acc_per_MPI;

    for ( esint i = 0; i < maxDevNumber; ++i ) {
#pragma omp parallel for 
        for ( esint d = 0 ; d < accDomains[i].size(); ++d ) {
            esint domN = accDomains[i].at(d);
            for ( esint j = 0 ; j < domains[domN]->K.cols; ++j ) {
                SparseKPack[i].SetX(d, j, (*tm1[domN]).at(j));   
            }
        }
    }

    int maxThreads = omp_get_max_threads();
    bool resetNested = false;

    if ( omp_get_max_active_levels() == 1 ) {
        omp_set_nested(1);
        resetNested = true;
    }


#pragma omp parallel num_threads( maxDevNumber + 1 ) 
    {
        int thread = omp_get_thread_num();
        if ( thread < maxDevNumber ) {
            MICtime[ thread ] = Measure::time();
            SparseKPack[ thread ].SolveMIC();
            esint end = (esint) accDomains[ thread ].size() * 
                SparseKPack[thread].getMICratio();
            for (esint i = 0 ; i < end; ++i) {
                esint domN = accDomains[ thread ].at( i );
                SparseKPack[thread].GetY( i, *tm2[domN] );
            }
            MICtime[ thread ] = Measure::time() - MICtime[ thread ];
        } else {
            omp_set_num_threads( maxThreads - maxDevNumber );     

            double startCPU = Measure::time();

            for ( esint i = 0; i < maxDevNumber; ++i ) {
                esint start = (esint) (SparseKPack[ i ].getNMatrices() * 
                        SparseKPack[ i ].getMICratio());
#pragma omp parallel for
                for (esint d = start; d < SparseKPack[ i ].getNMatrices(); ++d ) {
                    esint domN = accDomains[ i ].at( d );
                    domains[ domN ]->multKplusLocal(*tm1[domN], *tm2[ domN ]);
                }
            }
            
            for ( esint c = 0; c < numClusters; ++c ) {
                #pragma omp parallel for
                for ( esint d = 0; d < clusters[c].domains.size(); ++d ) {
                    esint e0_start = d*clusters[c].domains[d].Kplus_R.cols;
                    esint domain_size = clusters[c].domains[d].domain_prim_size;

                    clusters[c].domains[d].Kplus_R.DenseMatVec( clusters[c].vec_alfa, 
                        *tm3[clusters[c].domains[d].domain_global_index], 'N', e0_start);
                }
            }


//#pragma omp parallel for 
//            for (size_t d = 0; d < domains.size(); d++) {
//                esint e0_start	=  d	* domains[d]->Kplus_R.cols;
//                esint domain_size = domains[d]->domain_prim_size;
//
//                domains[d]->Kplus_R.DenseMatVec(vec_alfa, *tm3[d],'N', e0_start);
//            }

            CPUtime = Measure::time() - startCPU;
        }
    }

    if ( resetNested ) {
        omp_set_nested( 0 );
    }
    omp_set_num_threads( maxThreads );

#pragma omp parallel for
    for (size_t d = 0; d < domains.size(); d++) {
        esint domain_size = domains[d]->domain_prim_size;

        for (esint i = 0; i < domain_size; i++)
        {
            (*x_in[d])[i] = (*tm2[d])[i] + (*tm3[d])[i];
        }

    }

}


void SuperClusterAcc::init() {
    tm1.resize( number_of_subdomains_per_supercluster );
    tm2.resize( number_of_subdomains_per_supercluster );
    tm3.resize( number_of_subdomains_per_supercluster );

    for (esint c = 0; c < numClusters; c++) {
        // Get an original mapping of the subdomains
        for (int d = 0; d < clusters[c].domains.size(); d++) {
            tm1[clusters[c].domains[d].domain_global_index] = & clusters[c].tm1[d];
            tm2[clusters[c].domains[d].domain_global_index] = & clusters[c].tm2[d];
            tm3[clusters[c].domains[d].domain_global_index] = & clusters[c].tm3[d];
        }

    }
}


