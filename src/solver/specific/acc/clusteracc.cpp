
#include "clusteracc.h"

using namespace espreso;

ClusterAcc::~ClusterAcc() {
    if (this->deleteMatrices) {
        for (eslocal i = 0; i < config::solver::N_MICS; i++) {
            if (matricesPerAcc[i]) {
                //delete [] matricesPerAcc[i];
            }
        }
    }
}

void ClusterAcc::Create_SC_perDomain(bool USE_FLOAT) {

    ESINFO(PROGRESS2) << "Creating B1*K+*B1t : using MKL Pardiso on Xeon Phi accelerator : ";

    // First, get the available memory on coprocessors (in bytes)
    double usableRAM = 0.9;
    long micMem[config::solver::N_MICS];
    for (eslocal i = 0; i < config::solver::N_MICS; ++i) {
        long currentMem = 0;
#pragma offload target(mic:i)
        {
            long pages = sysconf(_SC_AVPHYS_PAGES);
            long page_size = sysconf(_SC_PAGE_SIZE);
            currentMem = pages * page_size;
        }
        micMem[i] = currentMem;
    }


    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);
    }

    // compute sizes of data to be offloaded to MIC
    eslocal maxDevNumber = config::solver::N_MICS;
    eslocal matrixPerPack[config::solver::N_MICS];
    for (eslocal i = 0 ; i < config::solver::N_MICS; ++i) {
        matrixPerPack[i] = domains.size() / maxDevNumber;
    }
    for (eslocal i = 0 ; i < domains.size() % maxDevNumber; ++i) {
        matrixPerPack[i]++;
    }
    eslocal offset = 0;
    bool symmetric = true;
    this->B1KplusPacks.resize( maxDevNumber );
    this->accDomains.resize( maxDevNumber );

    for ( int i = 0; i < maxDevNumber; i++ )
    {
        long dataSize = 0;
        long currDataSize = 0;
        bool MICFull = false;
        eslocal numCPUDomains = 0;

        for ( int j = offset; j < offset + matrixPerPack[i]; j++ ) {
            if (!symmetric) {
                currDataSize = domains[j].B1t_comp_dom.cols * domains[j].B1t_comp_dom.cols;
            } else {
                // isPacked => is symmetric
                currDataSize = ( ( 1.0 + ( double ) domains[j].B1t_comp_dom.cols ) *
                        ( ( double ) domains[j].B1t_comp_dom.cols ) / 2.0 );
            }
            long dataInBytes = (currDataSize + dataSize) * sizeof(double);
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

        this->B1KplusPacks[i].Resize( matrixPerPack[i], dataSize );

        for ( eslocal j = 0; j < matrixPerPack[i]; ++j ) {
            this->B1KplusPacks[ i ].PreparePack( j, domains[accDomains[i].at(j)].B1t_comp_dom.cols,
                    domains[accDomains[i].at(j)].B1t_comp_dom.cols, symmetric );
        }
    }

    bool assembleOnCPU = true;
    if ( assembleOnCPU ) {
        // iterate over available MICs and assemble their domains on CPU
        for (eslocal i = 0 ; i < config::solver::N_MICS ; ++i) {
            this->B1KplusPacks[i].AllocateVectors( );
            this->B1KplusPacks[i].SetDevice( i );
#pragma omp parallel for
            for (eslocal j = 0 ; j < accDomains[i].size() ; ++j) {
                eslocal domN = accDomains[i].at(j);
                double * matrixPointer = this->B1KplusPacks[i].getMatrixPointer(j);
                SparseMatrix tmpB;
                domains[domN].B1_comp_dom.MatTranspose(tmpB);
                SparseSolverCPU tmpsps;
                tmpsps.Create_SC_w_Mat( domains[domN].K, tmpB, domains[domN].B1Kplus, false, symmetric);
                memcpy(matrixPointer, &(domains[domN].B1Kplus.dense_values[0]), 
                        this->B1KplusPacks[i].getDataLength(j) * sizeof(double) );
                SEQ_VECTOR<double>().swap(  domains[domN].B1Kplus.dense_values);
                //domains[domN].B1Kplus.Clear();
                std::cout << ".";
            }
        }

    #pragma omp parallel num_threads(config::solver::N_MICS)
    {
         this->B1KplusPacks[omp_get_thread_num()].CopyToMIC();
    }

#pragma omp parallel for
        for (eslocal d = 0; d < hostDomains.size(); ++d) {

            std::cout << "*";
            eslocal domN = hostDomains.at(d);
            SparseMatrix TmpB;
            domains[domN].B1_comp_dom.MatTranspose(TmpB);
            SparseSolverCPU tmpsps;
            tmpsps.Create_SC_w_Mat( domains[domN].K, TmpB, domains[domN].B1Kplus, false, symmetric);
        }

    } else {
#pragma omp parallel     
        {
            if ((accDomains[omp_get_thread_num()].size()>0) && (omp_get_thread_num() < config::solver::N_MICS)) {
                // the first config::solver::N_MICS threads will communicate with MICs
                // now copy data to MIC and assemble Schur complements
                int  i = omp_get_thread_num();
                this->B1KplusPacks[i].AllocateVectors( );
                this->B1KplusPacks[i].SetDevice( i );
                SparseSolverAcc tmpsps_mic;
                SparseMatrix** K = new SparseMatrix*[matrixPerPack[i]];
                SparseMatrix** B = new SparseMatrix*[matrixPerPack[i]];

                for (eslocal j = 0; j < accDomains[i].size(); ++j) {
                    K[j] = &(domains[accDomains[i].at(j)].K);
                    B[j] = &(domains[accDomains[i].at(j)].B1t_comp_dom);
                }
                this->B1KplusPacks[i].CopyToMIC();
                tmpsps_mic.Create_SC_w_Mat( K, B, this->B1KplusPacks[i], accDomains[i].size(),  0, i);
            } else {
                // the remaining threads will assemble SC on CPU if necessary
#pragma omp single nowait
                {
                    for (eslocal d = 0; d < hostDomains.size(); ++d) {
#pragma omp task
                        {
                            //std::cout << ".";
                            eslocal domN = hostDomains.at(d);
                            SparseMatrix TmpB;
                            domains[domN].B1_comp_dom.MatTranspose(TmpB);
                            SparseSolverCPU tmpsps;
                            tmpsps.Create_SC_w_Mat( domains[domN].K, TmpB, domains[domN].B1Kplus, false, 0);
                        }
                    }
                }
            }
        }
    }


    /*
       cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
       domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

       if (cluster_global_index == 1)
       cout << "Creating B1*K+*B1t : using MKL Pardiso on Xeon Phi accelerator : ";
       this->NUM_MICS = 2;

    // compute sizes of data to be offloaded to MIC
    eslocal maxDevNumber = this->NUM_MICS;
    if (this->NUM_MICS == 0) {
    maxDevNumber = 1;
    }
    eslocal matrixPerPack = domains.size() / maxDevNumber;
    eslocal offset = 0;
    bool symmetric = true;
    this->B1KplusPacks.resize( maxDevNumber );
    eslocal * dom2dev = new eslocal[ domains.size() ];
    eslocal * offsets = new eslocal[maxDevNumber];

    for ( eslocal i = 0; i < maxDevNumber; i++ ) {
    if ( i == maxDevNumber - 1 ) {
    matrixPerPack += domains.size() % maxDevNumber;
    }

    long dataSize = 0;
    offsets[i] = offset;

    for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
    if (!symmetric) {
    dataSize += domains[j].B1t_comp_dom.cols * domains[j].B1t_comp_dom.cols;
    } else {
    // isPacked => is symmetric
    dataSize += ( ( 1.0 + ( double ) domains[j].B1t_comp_dom.cols ) *
    ( ( double ) domains[j].B1t_comp_dom.cols ) / 2.0 );
    }
    dom2dev[ j ] = i;
    }

    this->B1KplusPacks[i].Resize( matrixPerPack, dataSize );

    for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
    this->B1KplusPacks[ i ].PreparePack( j - offset, domains[j].B1t_comp_dom.cols,
    domains[j].B1t_comp_dom.cols,  symmetric );
    }
    offset += matrixPerPack;
    }

    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

    if (cluster_global_index == 1) cout << "."; // << i ;

    SparseSolverCPU tmpsps;
    if ( i == 0 && cluster_global_index == 1) tmpsps.msglvl = 1;
    tmpsps.Create_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B1Kplus, false, 1 );

    if (USE_FLOAT){
    domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
    domains[i].B1Kplus.USE_FLOAT = true;
    }

    this->B1KplusPacks[ dom2dev[ i ] ].AddDenseMatrix( i - offsets[dom2dev[i]], &(domains[i].B1Kplus.dense_values[0]) );
    domains[i].B1Kplus.Clear();

    }

    delete [] dom2dev;
    delete [] offsets;
    if (this->NUM_MICS == 0) {
    this->B1KplusPacks[0].AllocateVectors( );
    }
    for (eslocal i = 0; i < this->NUM_MICS ; i++ ) {
        this->B1KplusPacks[i].AllocateVectors( );
        this->B1KplusPacks[i].SetDevice( i );
        this->B1KplusPacks[i].CopyToMIC();
    }

    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1t_comp_dom.Clear();

    if (cluster_global_index == 1)
        cout << endl;
    */
}

void ClusterAcc::Create_Kinv_perDomain() {

    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

    ESINFO(PROGRESS2) << "Creating B1*K+*B1t on Xeon Phi accelerator";


    // compute sizes of data to be offloaded to MIC
    eslocal maxDevNumber = config::solver::N_MICS;
    if (config::solver::N_MICS == 0) {
        maxDevNumber = 1;
    }
    eslocal matrixPerPack = domains.size() / maxDevNumber;
    eslocal offset = 0;
    bool symmetric = true;
    this->B1KplusPacks.resize( maxDevNumber );
    eslocal * dom2dev = new eslocal[ domains.size() ];
    eslocal * offsets = new eslocal[maxDevNumber];

    for ( eslocal i = 0; i < maxDevNumber; i++ ) {
        if ( i == maxDevNumber - 1 ) {
            matrixPerPack += domains.size() % maxDevNumber;
        }

        long dataSize = 0;
        offsets[i] = offset;

        for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
            if (!symmetric) {
                dataSize += domains[j].B1t_comp_dom.cols * domains[j].B1t_comp_dom.cols;
            } else {
                // isPacked => is symmetric
                dataSize += ( ( 1.0 + ( double ) domains[j].B1t_comp_dom.cols ) *
                        ( ( double ) domains[j].B1t_comp_dom.cols ) / 2.0 );
            }
            dom2dev[ j ] = i;
        }

        this->B1KplusPacks[i].Resize( matrixPerPack, dataSize );
        for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
            this->B1KplusPacks[ i ].PreparePack( j - offset, domains[j].B1t_comp_dom.cols,
                    domains[j].B1t_comp_dom.cols,  symmetric );

        }
        offset += matrixPerPack;
    }

    ESINFO(PROGRESS2) << "Creating B1*K+*B1t : ";

    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

        domains[i].KplusF.msglvl = 0;

        if ( i == 0 && cluster_global_index == 1) {
            domains[i].KplusF.msglvl = Info::report(LIBRARIES) ? 1 : 0;
        }

        //SolveMatF is obsolete - use Schur complement instead
        domains[i].KplusF.SolveMatF(domains[i].B1t_comp_dom, domains[i].B1Kplus, false);
        domains[i].B1Kplus.MatTranspose();

        ESINFO(PROGRESS2) << Info::plain() << " " << i;

        SparseMatrix Btmp;
        Btmp.MatAddInPlace(domains[i].B1Kplus, 'N', 1.0);

        domains[i].B1Kplus.Clear ();
        domains[i].B1Kplus.MatMat(Btmp,'N', domains[i].B1t_comp_dom);
        domains[i].B1Kplus.ConvertCSRToDense(0);
        //domains[i].B1Kplus.ConvertDenseToDenseFloat(0);


        this->B1KplusPacks[ dom2dev[ i ] ].AddDenseMatrix( i - offsets[dom2dev[i]], &(domains[i].B1Kplus.dense_values[0]) );
        domains[i].B1Kplus.Clear();

    }

    ESINFO(PROGRESS2);

    delete [] dom2dev;
    delete [] offsets;
    if (config::solver::N_MICS == 0) {
        this->B1KplusPacks[0].AllocateVectors( );
    }
    for (eslocal i = 0; i < config::solver::N_MICS ; i++ ) {
        this->B1KplusPacks[i].AllocateVectors( );
        this->B1KplusPacks[i].SetDevice( i );
        this->B1KplusPacks[i].CopyToMIC();
    }


    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1t_comp_dom.Clear();

    //	std::cout << "This function is obsolete - use Create_SC_perDomain" << std::endl;
    //	return();

}

//void ClusterAcc::Create_SC_perDomain(bool USE_FLOAT) {
//
//    if (cluster_global_index == 1)
//        cout << "Creating B1*K+*B1t : using MKL Pardiso on Xeon Phi accelerator : ";
//
//    // First, get the available memory on coprocessors (in bytes)
//    long micMem[config::solver::N_MICS];
//    for (eslocal i = 0; i < config::solver::N_MICS; ++i) {
//        long currentMem = 0;
//        #pragma offload target(mic:0)
//        {
//            long pages = sysconf(_SC_AVPHYS_PAGES);
//            long page_size = sysconf(_SC_PAGE_SIZE);
//            currentMem = pages * page_size;
//        }
//        micMem[i] = currentMem;
//    }
//    for (eslocal i = 0; i < config::solver::N_MICS; ++i) {
//        totalMem += micMem[i];
//    }
//
//
//    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {
//        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);
//    }
//
//    // compute sizes of data to be offloaded to MIC
//    int maxDevNumber = config::solver::N_MICS;
//    if (config::solver::N_MICS == 0) {
//        maxDevNumber = 1;
//    }
//    int matrixPerPack = domains.size() / maxDevNumber;
//    int offset = 0;
//    bool symmetric = true;
//    this->B1KplusPacks.resize( maxDevNumber );
//    int * dom2dev = new int[ domains.size() ];
//    int * offsets = new int[ maxDevNumber ];
//
//    for ( int i = 0; i < maxDevNumber; i++ )
//    {
//        if ( i == maxDevNumber - 1 ) {
//            matrixPerPack += domains.size() % maxDevNumber;
//        }
//
//        long dataSize = 0;
//        offsets[i] = offset;
//
//        for ( int j = offset; j < offset + matrixPerPack; j++ ) {
//            if (!symmetric) {
//                dataSize += domains[j].B1t_comp_dom.cols * domains[j].B1t_comp_dom.cols;
//            } else {
//                // isPacked => is symmetric
//                dataSize += ( ( 1.0 + ( double ) domains[j].B1t_comp_dom.cols ) *
//                        ( ( double ) domains[j].B1t_comp_dom.cols ) / 2.0 );
//            }
//            dom2dev[ j ] = i;
//        }
//
//        this->B1KplusPacks[i].Resize( matrixPerPack, dataSize );
//
//        for ( int j = offset; j < offset + matrixPerPack; j++ ) {
//            this->B1KplusPacks[ i ].PreparePack( j - offset, domains[j].B1t_comp_dom.cols,
//                    domains[j].B1t_comp_dom.cols,  symmetric );
//        }
//        offset += matrixPerPack;
//    }
//    int matrixPP = domains.size() / maxDevNumber;
//#pragma omp parallel num_threads(config::solver::N_MICS)
//    {
//        int  i = omp_get_thread_num();
//        std::cout << "DEVICE: " <<i << std::endl;
//        this->B1KplusPacks[i].AllocateVectors( );
//        this->B1KplusPacks[i].SetDevice( i );
//        SparseSolverAcc tmpsps_mic;
//        if ( i == maxDevNumber - 1 ) {
//            matrixPP += domains.size() % maxDevNumber;
//        }
//        SparseMatrix** K = new SparseMatrix*[matrixPP];
//        SparseMatrix** B = new SparseMatrix*[matrixPP];
//
//        for (int j = offsets[i]; j < offsets[i] + matrixPP; j++ ) {
//            K[j-offsets[i]] = &(domains[j].K);
//            B[j-offsets[i]] = &(domains[j].B1t_comp_dom);
//        }
//        this->B1KplusPacks[i].CopyToMIC();
//        tmpsps_mic.Create_SC_w_Mat( K, B, this->B1KplusPacks[i], matrixPP,  0, i);
//      //  MKL_INT * SC_sizes = new MKL_INT[matrixPP];
//      //  for (eslocal j  = 0; j < matrixPP; j++) {
//      //      SC_sizes[j] = domains[offsets[i]+j].B1t_comp_dom.cols;
//      //  }
//      //  this->solver[i].Create_SC(this->B1KplusPacks[i],SC_sizes, 0);
//    }
//    delete [] dom2dev;
//    delete [] offsets;
//    /*
//       cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
//       domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);
//
//       if (cluster_global_index == 1)
//       cout << "Creating B1*K+*B1t : using MKL Pardiso on Xeon Phi accelerator : ";
//       this->NUM_MICS = 2;
//
//    // compute sizes of data to be offloaded to MIC
//    eslocal maxDevNumber = this->NUM_MICS;
//    if (this->NUM_MICS == 0) {
//    maxDevNumber = 1;
//    }
//    eslocal matrixPerPack = domains.size() / maxDevNumber;
//    eslocal offset = 0;
//    bool symmetric = true;
//    this->B1KplusPacks.resize( maxDevNumber );
//    eslocal * dom2dev = new eslocal[ domains.size() ];
//    eslocal * offsets = new eslocal[maxDevNumber];
//
//    for ( eslocal i = 0; i < maxDevNumber; i++ ) {
//    if ( i == maxDevNumber - 1 ) {
//    matrixPerPack += domains.size() % maxDevNumber;
//    }
//
//    long dataSize = 0;
//    offsets[i] = offset;
//
//    for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
//    if (!symmetric) {
//    dataSize += domains[j].B1t_comp_dom.cols * domains[j].B1t_comp_dom.cols;
//    } else {
//    // isPacked => is symmetric
//    dataSize += ( ( 1.0 + ( double ) domains[j].B1t_comp_dom.cols ) *
//    ( ( double ) domains[j].B1t_comp_dom.cols ) / 2.0 );
//    }
//    dom2dev[ j ] = i;
//    }
//
//    this->B1KplusPacks[i].Resize( matrixPerPack, dataSize );
//
//    for ( eslocal j = offset; j < offset + matrixPerPack; j++ ) {
//    this->B1KplusPacks[ i ].PreparePack( j - offset, domains[j].B1t_comp_dom.cols,
//    domains[j].B1t_comp_dom.cols,  symmetric );
//    }
//    offset += matrixPerPack;
//    }
//
//    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {
//
//    if (cluster_global_index == 1) cout << "."; // << i ;
//
//    SparseSolverCPU tmpsps;
//    if ( i == 0 && cluster_global_index == 1) tmpsps.msglvl = 1;
//    tmpsps.Create_SC_w_Mat( domains[i].K, domains[i].B1t_comp_dom, domains[i].B1Kplus, false, 1 );
//
//    if (USE_FLOAT){
//    domains[i].B1Kplus.ConvertDenseToDenseFloat( 1 );
//    domains[i].B1Kplus.USE_FLOAT = true;
//    }
//
//    this->B1KplusPacks[ dom2dev[ i ] ].AddDenseMatrix( i - offsets[dom2dev[i]], &(domains[i].B1Kplus.dense_values[0]) );
//    domains[i].B1Kplus.Clear();
//
//    }
//
//    delete [] dom2dev;
//    delete [] offsets;
//    if (this->NUM_MICS == 0) {
//    this->B1KplusPacks[0].AllocateVectors( );
//    }
//    for (eslocal i = 0; i < this->NUM_MICS ; i++ ) {
//        this->B1KplusPacks[i].AllocateVectors( );
//        this->B1KplusPacks[i].SetDevice( i );
//        this->B1KplusPacks[i].CopyToMIC();
//    }
//
//    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
//        domains[i].B1t_comp_dom.Clear();
//
//    if (cluster_global_index == 1)
//        cout << endl;
//    */
//#pragma offload target(mic:0)
//    {
//        long pages = sysconf(_SC_AVPHYS_PAGES);
//        long page_size = sysconf(_SC_PAGE_SIZE);
//        std::cout << pages * page_size << std::endl;
//    }
//}


//TODO:
void ClusterAcc::SetupKsolvers ( ) {
    // this part is for setting CPU pardiso, temporarily until everything is
    // solved on MIC
    cilk_for (eslocal d = 0; d < domains.size(); d++) {

        // Import of Regularized matrix K into Kplus (Sparse Solver)
        switch (config::solver::KSOLVER) {
            case 0: {
                        domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
                        break;
                    }
            case 1: {
                        domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
                        break;
                    }
            case 2: {
                        domains[d].Kplus.ImportMatrix_fl(domains[d].K);
                        break;
                    }
            case 3: {
                        domains[d].Kplus.ImportMatrix_fl(domains[d].K);
                        break;
                    }
            case 4: {
                        domains[d].Kplus.ImportMatrix_fl(domains[d].K);
                        break;
                    }
            default:
                    ESINFO(ERROR) << "Invalid KSOLVER value.";
                    exit(EXIT_FAILURE);
        }

        if (config::solver::KEEP_FACTORS) {
            std::stringstream ss;
            ss << "init -> rank: " << config::env::MPIrank << ", subdomain: " << d;
            domains[d].Kplus.keep_factors = true;
            if (config::solver::KSOLVER != 1) {
                domains[d].Kplus.Factorization (ss.str());
            }
        } else {
            domains[d].Kplus.keep_factors = false;
            domains[d].Kplus.MPIrank = config::env::MPIrank;
        }

        domains[d].domain_prim_size = domains[d].Kplus.cols;

        if ( d == 0 && config::env::MPIrank == 0) {
            domains[d].Kplus.msglvl = Info::report(LIBRARIES) ? 1 : 0;
        }
    }
    if (!USE_KINV) {
        // send matrices to Xeon Phi
        eslocal nMatrices = domains.size();
        this->matricesPerAcc.reserve(config::solver::N_MICS);
        SEQ_VECTOR<eslocal> nMatPerMIC;
        nMatPerMIC.resize(config::solver::N_MICS);

        for (eslocal i = 0; i < config::solver::N_MICS; i++) {
            nMatPerMIC[i] = nMatrices / config::solver::N_MICS;
        }

        for (eslocal i = 0 ; i < nMatrices % config::solver::N_MICS; i++ ) {
            nMatPerMIC[i]++;
        }

        eslocal offset = 0;
        for (eslocal i = 0; i < config::solver::N_MICS; i++) {
            for (eslocal j = offset; j < offset + nMatPerMIC[ i ]; j++) {
                (this->matricesPerAcc[i])[j - offset] = &(domains[j].K);
            }
            offset += nMatPerMIC[i];
        }
        this->solver.resize(config::solver::N_MICS);

#pragma omp parallel num_threads(config::solver::N_MICS)
        {
            eslocal myAcc = omp_get_thread_num();
            this->solver[myAcc].ImportMatrices_wo_Copy(this->matricesPerAcc[myAcc], nMatPerMIC[myAcc], myAcc);
            this->solver[myAcc].Factorization("");
        }


        deleteMatrices = true  ;
    }
}
