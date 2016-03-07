
#include "clusteracc.h"

ClusterAcc::~ClusterAcc() {
    if (this->deleteMatrices) {
        for (eslocal i = 0; i < this->NUM_MICS; i++) {
            if (matricesPerAcc[i]) {
                //delete [] matricesPerAcc[i];
            }
        }
    }
}

void ClusterAcc::Create_SC_perDomain(bool USE_FLOAT) {

    if (cluster_global_index == 1)
        cout << "Creating B1*K+*B1t : using MKL Pardiso on Xeon Phi accelerator : ";

    this->NUM_MICS = 2;

    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);
    }

    // compute sizes of data to be offloaded to MIC
    int maxDevNumber = this->NUM_MICS;
    if (this->NUM_MICS == 0) {
        maxDevNumber = 1;
    }
    int matrixPerPack = domains.size() / maxDevNumber;
    int offset = 0;
    bool symmetric = true;
    this->B1KplusPacks.resize( maxDevNumber );
    int * dom2dev = new int[ domains.size() ];
    int * offsets = new int[ maxDevNumber ];

    for ( int i = 0; i < maxDevNumber; i++ )
    {
        if ( i == maxDevNumber - 1 ) {
            matrixPerPack += domains.size() % maxDevNumber;
        }

        long dataSize = 0;
        offsets[i] = offset;

        for ( int j = offset; j < offset + matrixPerPack; j++ ) {
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

        for ( int j = offset; j < offset + matrixPerPack; j++ ) {
            this->B1KplusPacks[ i ].PreparePack( j - offset, domains[j].B1t_comp_dom.cols,
                    domains[j].B1t_comp_dom.cols,  symmetric );
        }
        offset += matrixPerPack;
    }
    int matrixPP = domains.size() / maxDevNumber;
#pragma omp parallel num_threads(NUM_MICS)
    {
        int  i = omp_get_thread_num();
        std::cout << "DEVICE: " <<i << std::endl;
        this->B1KplusPacks[i].AllocateVectors( );
        this->B1KplusPacks[i].SetDevice( i );
        SparseSolverAcc tmpsps_mic;
        if ( i == maxDevNumber - 1 ) {
            matrixPP += domains.size() % maxDevNumber;
        }
        SparseMatrix** K = new SparseMatrix*[matrixPP];
        SparseMatrix** B = new SparseMatrix*[matrixPP];

        for (int j = offsets[i]; j < offsets[i] + matrixPP; j++ ) {
            K[j-offsets[i]] = &(domains[j].K);
            B[j-offsets[i]] = &(domains[j].B1t_comp_dom);
        }
        this->B1KplusPacks[i].CopyToMIC();
        tmpsps_mic.Create_SC_w_Mat( K, B, this->B1KplusPacks[i], matrixPP,  0, i);
    }

    delete [] dom2dev;
    delete [] offsets;
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

    if (cluster_global_index == 1)
        cout << "Creating B1*K+*B1t on Xeon Phi accelerator : ";

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

        domains[i].KplusF.msglvl = 0;

        if ( i == 0 && cluster_global_index == 1) domains[i].KplusF.msglvl=1;

        //SolveMatF is obsolete - use Schur complement instead
        domains[i].KplusF.SolveMatF(domains[i].B1t_comp_dom, domains[i].B1Kplus, false);
        domains[i].B1Kplus.MatTranspose();

        if (cluster_global_index == 1 && i == 0)
            cout << "Creating B1*K+*B1t : ";

        if (cluster_global_index == 1) {
            cout << " " << i ;
        }

        SparseMatrix Btmp;
        Btmp.MatAddInPlace(domains[i].B1Kplus, 'N', 1.0);

        domains[i].B1Kplus.Clear ();
        domains[i].B1Kplus.MatMat(Btmp,'N', domains[i].B1t_comp_dom);
        domains[i].B1Kplus.ConvertCSRToDense(0);
        //domains[i].B1Kplus.ConvertDenseToDenseFloat(0);


        this->B1KplusPacks[ dom2dev[ i ] ].AddDenseMatrix( i - offsets[dom2dev[i]], &(domains[i].B1Kplus.dense_values[0]) );
        domains[i].B1Kplus.Clear();

    }

    cout << endl;

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

    //	std::cout << "This function is obsolete - use Create_SC_perDomain" << std::endl;
    //	return();

}



//TODO:
void ClusterAcc::SetupKsolvers ( ) {
    // this part is for setting CPU pardiso, temporarily until everything is
    // solved on MIC
    cilk_for (eslocal d = 0; d < domains.size(); d++) {

        // Import of Regularized matrix K into Kplus (Sparse Solver)
        switch (esconfig::solver::KSOLVER) {
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
                    ESLOG(eslog::ERROR) << "Invalid KSOLVER value.";
                    exit(EXIT_FAILURE);
        }

        if (esconfig::solver::KEEP_FACTORS == 1) {
            std::stringstream ss;
            ss << "init -> rank: " << esconfig::MPIrank << ", subdomain: " << d;
            domains[d].Kplus.keep_factors = true;
            if (esconfig::solver::KSOLVER != 1) {
                domains[d].Kplus.Factorization (ss.str());
            }
        } else {
            domains[d].Kplus.keep_factors = false;
            domains[d].Kplus.MPIrank = esconfig::MPIrank;
        }

        domains[d].domain_prim_size = domains[d].Kplus.cols;

        if ( d == 0 && esconfig::MPIrank == 0) domains[d].Kplus.msglvl=0;
        if (esconfig::MPIrank == 0) std::cout << ".";

    }
std::cout << "TEST 0" <<std::endl;
    // send matrices to Xeon Phi
    eslocal nMatrices = domains.size();  
    this->matricesPerAcc.reserve(this->NUM_MICS);
    SEQ_VECTOR<eslocal> nMatPerMIC;
    nMatPerMIC.resize(NUM_MICS);
    
    for (eslocal i = 0; i < this->NUM_MICS; i++) {
        nMatPerMIC[i] = nMatrices / this->NUM_MICS;
    }
    
    for (eslocal i = 0 ; i < nMatrices % this->NUM_MICS; i++ ) {
        nMatPerMIC[i]++;
    }
    
    eslocal offset = 0;
    for (eslocal i = 0; i < this->NUM_MICS; i++) {
        this->matricesPerAcc[i] = new SparseMatrix*[ nMatPerMIC[ i ] ];
        for (eslocal j = offset; j < offset + nMatPerMIC[ i ]; j++) {
            (this->matricesPerAcc[i])[j - offset] = &(domains[j].K); 
        }
        offset += nMatPerMIC[i];
    }
    this->solver.resize(this->NUM_MICS);
    for (eslocal i = 0; i < this->NUM_MICS; i++) {
        this->solver[i].ImportMatrices_wo_Copy(this->matricesPerAcc[i], nMatPerMIC[i], i);
        this->solver[i].Factorization("");
    }

    

    deleteMatrices = true  ;
}
