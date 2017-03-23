
#include "clusteracc.h"

#include "../../../assembler/instance.h"

using namespace espreso;

ClusterAcc::~ClusterAcc() {
    if (this->deleteMatrices) {
        for (eslocal i = 0; i < configuration.N_MICS; i++) {
            if (matricesPerAcc[i]) {
                //delete [] matricesPerAcc[i];
            }
        }
    }
}

void ClusterAcc::Create_SC_perDomain(bool USE_FLOAT) {

    ESINFO(PROGRESS3) << "Creating Local Schur complements";

    // detect how many MPI processes is running per node
    int _MPIglobalRank;
    int _MPIglobalSize;
    int _MPInodeRank;
    int _MPInodeSize;
    std::string _nodeName;
    MPI_Comm _currentNode;
    MPI_Comm _storingProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &_MPIglobalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &_MPIglobalSize);
    int color;
    int length;
    std::vector<char> name(MPI_MAX_PROCESSOR_NAME);
    MPI_Get_processor_name(name.data(), &length);
    _nodeName = std::string(name.begin(), name.begin() + length);
    std::vector<int> rCounts(_MPIglobalSize);
    std::vector<int> rDispl(_MPIglobalSize);
    std::vector<int> colors(_MPIglobalSize);
    std::vector<char> names;
    MPI_Gather(&length, sizeof(int), MPI_BYTE, rCounts.data(), sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
    for (size_t i = 1; i < rCounts.size(); i++) {
        rDispl[i] += rDispl[i - 1] + rCounts[i - 1];
    }
    names.resize(rDispl.back() + rCounts.back());
    MPI_Gatherv(name.data(), length * sizeof(char), MPI_BYTE, names.data(), rCounts.data(), rDispl.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
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
    MPI_Scatter(colors.data(), sizeof(int), MPI_BYTE, &color, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Comm_split(MPI_COMM_WORLD, color, _MPIglobalRank, &_currentNode);
    MPI_Comm_rank(_currentNode, &_MPInodeRank);
    MPI_Comm_size(_currentNode, &_MPInodeSize);
    MPI_Comm_split(MPI_COMM_WORLD, _MPInodeRank, _MPIglobalRank, &_storingProcs);
    // END - detect how many MPI processes is running per node
    ESINFO(PROGRESS3) << "MPI ranks per node: " << _MPInodeSize;
    //TODO:
    int accelerators_per_node = 2;
    int target;
    
    if (_MPInodeRank < (_MPInodeSize/2)) //  .../accelerators_per_node
        target = 0;
    else
        target = 1;
    // Set process placement on Xeon Phi
    int ranks_per_node = _MPInodeSize;
    int rank = _MPInodeRank % (_MPInodeSize/2); //  .../accelerators_per_node
    for (eslocal i = 0; i < configuration.N_MICS; ++i) {
        int used_core_num = 0;
        int first_core = 0;
        int last_core = 0;
        // TODO: #pragma offload target(mic:i)
#pragma offload target(mic:target)
        {
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
            cpu_set_t my_set;        /* Define your cpu_set bit mask. */
            CPU_ZERO(&my_set);       /* Initialize it all to 0, i.e. no CPUs selected. */
            int cores = sysconf(_SC_NPROCESSORS_ONLN); // for Xeon Phi 7120 - results is 244
            cores = (cores / 4 ) - 1; // keep core for system and remove effect of hyperthreading
            //std::cout << sysconf(_SC_NPROCESSORS_ONLN) << " Cores\n";//std::thread::hardware_concurrency() << "Cores \n";
            int cores_per_rank = accelerators_per_node * cores/ranks_per_node;
            for (int i = 0; i < cores_per_rank; i++) {
                if (i == 0) {
                    first_core = 1*(cores_per_rank)*rank + 1*i;
                }
                last_core = 1*(cores_per_rank)*rank + 1*i;
                int core = 1 + 4*(cores_per_rank)*rank + 4*i;
                CPU_SET(core, &my_set);     /* set the bit that represents core 7. */
                used_core_num++;
            }
            omp_set_num_threads(cores_per_rank-0);
            sched_setaffinity(0, sizeof(cpu_set_t), &my_set); /* Set affinity of tihs process to */
            /* the defined mask, i.e. only 7. */
        }
        //ESINFO(PROGRESS3)
        std::cout << "Global MPI rank: " << _MPIglobalRank << " - Node MPI rank: " << _MPInodeRank << " uses: " << used_core_num << " Xeon Phi processing cores of accelerator #" << target <<" (from " << first_core << " to " << last_core <<  ")\n";
    }


    // Ratio of work done on MIC
    double MICr = 1.0;
    if ( configuration.load_balancing ) {
        MICr = 0.1;
    }

    // First, get the available memory on coprocessors (in bytes)
    double usableRAM = 0.9;
    long micMem[configuration.N_MICS];
    for (eslocal i = 0; i < configuration.N_MICS; ++i) {
        long currentMem = 0;
//#pragma offload target(mic:i)
#pragma offload target(mic:target)
        {
            long pages = sysconf(_SC_AVPHYS_PAGES);
            long page_size = sysconf(_SC_PAGE_SIZE);
            currentMem = pages * page_size;
        }
        micMem[i] = currentMem;
    }


    #pragma omp parallel for
for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);
    }

    // compute sizes of data to be offloaded to MIC
    eslocal maxDevNumber = configuration.N_MICS;
    eslocal matrixPerPack[configuration.N_MICS];
    for (eslocal i = 0 ; i < configuration.N_MICS; ++i) {
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
        this->B1KplusPacks[i].setMICratio( MICr );

        if ( configuration.load_balancing ) {
            this->B1KplusPacks[i].enableLoadBalancing();
        } else {
            this->B1KplusPacks[i].disableLoadBalancing();
        }

        for ( eslocal j = 0; j < matrixPerPack[i]; ++j ) {
            this->B1KplusPacks[ i ].PreparePack( j, domains[accDomains[i].at(j)].B1t_comp_dom.cols,
                    domains[accDomains[i].at(j)].B1t_comp_dom.cols, symmetric );
        }
    }

    bool assembleOnCPU = true;
    if ( assembleOnCPU ) {
        // iterate over available MICs and assemble their domains on CPU
        for (eslocal i = 0 ; i < configuration.N_MICS ; ++i) {
            this->B1KplusPacks[i].AllocateVectors( );
            this->B1KplusPacks[i].SetDevice( target );
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
                ESINFO(PROGRESS3) << Info::plain() << ".";
            }
        }

    #pragma omp parallel num_threads(configuration.N_MICS)
    {
         this->B1KplusPacks[omp_get_thread_num()].CopyToMIC();
    }

#pragma omp parallel for
        for (eslocal d = 0; d < hostDomains.size(); ++d) {

            ESINFO(PROGRESS3) << Info::plain() << "*";
            eslocal domN = hostDomains.at(d);
            SparseMatrix TmpB;
            domains[domN].B1_comp_dom.MatTranspose(TmpB);
            SparseSolverCPU tmpsps;
            tmpsps.Create_SC_w_Mat( domains[domN].K, TmpB, domains[domN].B1Kplus, false, symmetric);
        }

    } else {
#pragma omp parallel
        {
            if ((accDomains[omp_get_thread_num()].size()>0) && (omp_get_thread_num() < configuration.N_MICS)) {
                // the first configuration.N_MICS threads will communicate with MICs
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

    #pragma omp parallel for
	for  (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

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

    #pragma omp parallel for
	for  (eslocal i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1t_comp_dom.Clear();

    if (cluster_global_index == 1)
        cout << endl;
    */
}

void ClusterAcc::Create_Kinv_perDomain() {

    #pragma omp parallel for
for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);

    ESINFO(PROGRESS3) << "Creating B1*K+*B1t on Xeon Phi accelerator";


    // compute sizes of data to be offloaded to MIC
    eslocal maxDevNumber = configuration.N_MICS;
    if (configuration.N_MICS == 0) {
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

    ESINFO(PROGRESS3) << "Creating B1*K+*B1t : ";

    #pragma omp parallel for
for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {

        domains[i].KplusF.msglvl = 0;

        if ( i == 0 && cluster_global_index == 1) {
            domains[i].KplusF.msglvl = Info::report(LIBRARIES) ? 1 : 0;
        }

        //SolveMatF is obsolete - use Schur complement instead
        domains[i].KplusF.SolveMatF(domains[i].B1t_comp_dom, domains[i].B1Kplus, false);
        domains[i].B1Kplus.MatTranspose();

        ESINFO(PROGRESS3) << Info::plain() << " " << i;

        SparseMatrix Btmp;
        Btmp.MatAddInPlace(domains[i].B1Kplus, 'N', 1.0);

        domains[i].B1Kplus.Clear ();
        domains[i].B1Kplus.MatMat(Btmp,'N', domains[i].B1t_comp_dom);
        domains[i].B1Kplus.ConvertCSRToDense(0);
        //domains[i].B1Kplus.ConvertDenseToDenseFloat(0);


        this->B1KplusPacks[ dom2dev[ i ] ].AddDenseMatrix( i - offsets[dom2dev[i]], &(domains[i].B1Kplus.dense_values[0]) );
        domains[i].B1Kplus.Clear();

    }

    ESINFO(PROGRESS3);

    delete [] dom2dev;
    delete [] offsets;
    if (configuration.N_MICS == 0) {
        this->B1KplusPacks[0].AllocateVectors( );
    }
    for (eslocal i = 0; i < configuration.N_MICS ; i++ ) {
        this->B1KplusPacks[i].AllocateVectors( );
        this->B1KplusPacks[i].SetDevice( i );
        this->B1KplusPacks[i].CopyToMIC();
    }


    #pragma omp parallel for
for (eslocal i = 0; i < domains_in_global_index.size(); i++ )
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
//    long micMem[configuration.N_MICS];
//    for (eslocal i = 0; i < configuration.N_MICS; ++i) {
//        long currentMem = 0;
//        #pragma offload target(mic:0)
//        {
//            long pages = sysconf(_SC_AVPHYS_PAGES);
//            long page_size = sysconf(_SC_PAGE_SIZE);
//            currentMem = pages * page_size;
//        }
//        micMem[i] = currentMem;
//    }
//    for (eslocal i = 0; i < configuration.N_MICS; ++i) {
//        totalMem += micMem[i];
//    }
//
//
//    cilk_for (eslocal i = 0; i < domains_in_global_index.size(); i++ ) {
//        domains[i].B1_comp_dom.MatTranspose(domains[i].B1t_comp_dom);
//    }
//
//    // compute sizes of data to be offloaded to MIC
//    int maxDevNumber = configuration.N_MICS;
//    if (configuration.N_MICS == 0) {
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
//#pragma omp parallel num_threads(configuration.N_MICS)
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
	#pragma omp parallel for
	for (eslocal d = 0; d < domains.size(); d++) {

        // Import of Regularized matrix K into Kplus (Sparse Solver)
		switch (configuration.Ksolver) {
		case ESPRESO_KSOLVER::DIRECT_DP:
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		case ESPRESO_KSOLVER::ITERATIVE:
			domains[d].Kplus.ImportMatrix_wo_Copy (domains[d].K);
			break;
		case ESPRESO_KSOLVER::DIRECT_SP:
			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
			break;
		case ESPRESO_KSOLVER::DIRECT_MP:
			domains[d].Kplus.ImportMatrix_fl(domains[d].K);
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
            ss << "init -> rank: " << environment->MPIrank << ", subdomain: " << d;
            domains[d].Kplus.keep_factors = true;
            if (configuration.Ksolver != ESPRESO_KSOLVER::ITERATIVE) {
                domains[d].Kplus.Factorization (ss.str());
            }
        } else {
            domains[d].Kplus.keep_factors = false;
            domains[d].Kplus.MPIrank = environment->MPIrank;
        }

        if ( d == 0 && environment->MPIrank == 0) {
            domains[d].Kplus.msglvl = Info::report(LIBRARIES) ? 1 : 0;
        }
    }
    if (!USE_KINV) {
        // send matrices to Xeon Phi
        eslocal nMatrices = domains.size();
        this->matricesPerAcc.reserve(configuration.N_MICS);
        SEQ_VECTOR<eslocal> nMatPerMIC;
        nMatPerMIC.resize(configuration.N_MICS);

        for (eslocal i = 0; i < configuration.N_MICS; i++) {
            nMatPerMIC[i] = nMatrices / configuration.N_MICS;
        }

        for (eslocal i = 0 ; i < nMatrices % configuration.N_MICS; i++ ) {
            nMatPerMIC[i]++;
        }

        eslocal offset = 0;
        for (eslocal i = 0; i < configuration.N_MICS; i++) {
            for (eslocal j = offset; j < offset + nMatPerMIC[ i ]; j++) {
                (this->matricesPerAcc[i])[j - offset] = &(domains[j].K);
            }
            offset += nMatPerMIC[i];
        }
        this->solver.resize(configuration.N_MICS);

#pragma omp parallel num_threads(configuration.N_MICS)
        {
            eslocal myAcc = omp_get_thread_num();
            this->solver[myAcc].ImportMatrices_wo_Copy(this->matricesPerAcc[myAcc], nMatPerMIC[myAcc], myAcc);
            this->solver[myAcc].Factorization("");
        }


        deleteMatrices = true  ;
    }
}

void ClusterAcc::CreateDirichletPrec( Instance *instance ) {
    // Prepare matrix pack

    // detect how many MPI processes is running per node
    int _MPIglobalRank;
    int _MPIglobalSize;
    int _MPInodeRank;
    int _MPInodeSize;
    std::string _nodeName;
    MPI_Comm _currentNode;
    MPI_Comm _storingProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &_MPIglobalRank);
    MPI_Comm_size(MPI_COMM_WORLD, &_MPIglobalSize);
    int color;
    int length;
    std::vector<char> name(MPI_MAX_PROCESSOR_NAME);
    MPI_Get_processor_name(name.data(), &length);
    _nodeName = std::string(name.begin(), name.begin() + length);
    std::vector<int> rCounts(_MPIglobalSize);
    std::vector<int> rDispl(_MPIglobalSize);
    std::vector<int> colors(_MPIglobalSize);
    std::vector<char> names;
    MPI_Gather(&length, sizeof(int), MPI_BYTE, rCounts.data(), sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
    for (size_t i = 1; i < rCounts.size(); i++) {
        rDispl[i] += rDispl[i - 1] + rCounts[i - 1];
    }
    names.resize(rDispl.back() + rCounts.back());
    MPI_Gatherv(name.data(), length * sizeof(char), MPI_BYTE, names.data(), rCounts.data(), rDispl.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
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
    MPI_Scatter(colors.data(), sizeof(int), MPI_BYTE, &color, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Comm_split(MPI_COMM_WORLD, color, _MPIglobalRank, &_currentNode);
    MPI_Comm_rank(_currentNode, &_MPInodeRank);
    MPI_Comm_size(_currentNode, &_MPInodeSize);
    MPI_Comm_split(MPI_COMM_WORLD, _MPInodeRank, _MPIglobalRank, &_storingProcs);
    // END - detect how many MPI processes is running per node
    ESINFO(PROGRESS3) << "MPI ranks per node: " << _MPInodeSize;
    //TODO:
    int accelerators_per_node = 2;
    int target;
    
    if (_MPInodeRank < (_MPInodeSize/2)) //  .../accelerators_per_node
        target = 0;
    else
        target = 1;
    // Set process placement on Xeon Phi
    int ranks_per_node = _MPInodeSize;
    int rank = _MPInodeRank % (_MPInodeSize/2); //  .../accelerators_per_node
    for (eslocal i = 0; i < configuration.N_MICS; ++i) {
        int used_core_num = 0;
        int first_core = 0;
        int last_core = 0;
        // TODO: #pragma offload target(mic:i)
#pragma offload target(mic:target)
        {
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
            cpu_set_t my_set;        /* Define your cpu_set bit mask. */
            CPU_ZERO(&my_set);       /* Initialize it all to 0, i.e. no CPUs selected. */
            int cores = sysconf(_SC_NPROCESSORS_ONLN); // for Xeon Phi 7120 - results is 244
            cores = (cores / 4 ) - 1; // keep core for system and remove effect of hyperthreading
            //std::cout << sysconf(_SC_NPROCESSORS_ONLN) << " Cores\n";//std::thread::hardware_concurrency() << "Cores \n";
            int cores_per_rank = accelerators_per_node * cores/ranks_per_node;
            for (int i = 0; i < cores_per_rank; i++) {
                if (i == 0) {
                    first_core = 1*(cores_per_rank)*rank + 1*i;
                }
                last_core = 1*(cores_per_rank)*rank + 1*i;
                int core = 1 + 4*(cores_per_rank)*rank + 4*i;
                CPU_SET(core, &my_set);     /* set the bit that represents core 7. */
                used_core_num++;
            }
            omp_set_num_threads(cores_per_rank-0);
            sched_setaffinity(0, sizeof(cpu_set_t), &my_set); /* Set affinity of tihs process to */
            /* the defined mask, i.e. only 7. */
        }
        //ESINFO(PROGRESS3)
        std::cout << "Global MPI rank: " << _MPIglobalRank << " - Node MPI rank: " << _MPInodeRank << " uses: " << used_core_num << " Xeon Phi processing cores of accelerator #" << target <<" (from " << first_core << " to " << last_core <<  ")\n";
    }

    // Ratio of work done on MIC
    double MICr = 1.0;
    if ( configuration.load_balancing ) {
        MICr = 0.1;
    }

    // First, get the available memory on coprocessors (in bytes)
    double usableRAM = 0.9;
    long micMem[configuration.N_MICS];
    for (eslocal i = 0; i < configuration.N_MICS; ++i) {
        long currentMem = 0;
//#pragma offload target(mic:i)
#pragma offload target(mic:target)
        {
            long pages = sysconf(_SC_AVPHYS_PAGES);
            long page_size = sysconf(_SC_PAGE_SIZE);
            currentMem = pages * page_size;
        }
        micMem[i] = currentMem;
    }

    // compute sizes of data to be offloaded to MIC
    eslocal maxDevNumber = configuration.N_MICS;
    eslocal matrixPerPack[configuration.N_MICS];
    for (eslocal i = 0 ; i < configuration.N_MICS; ++i) {
        matrixPerPack[i] = domains.size() / maxDevNumber;
    }
    for (eslocal i = 0 ; i < domains.size() % maxDevNumber; ++i) {
        matrixPerPack[i]++;
    }
    eslocal offset = 0;
    bool symmetric = true;
    this->DirichletPacks.resize( maxDevNumber );
    this->accPreconditioners.resize( maxDevNumber );

    for ( int i = 0; i < maxDevNumber; i++ )
    {
        long dataSize = 0;
        long currDataSize = 0;
        bool MICFull = false;
        eslocal numCPUDomains = 0;

        for ( int j = offset; j < offset + matrixPerPack[i]; j++ ) {
            eslocal size = domains[j].B1t_Dir_perm_vec.size();
            if (!symmetric) {
                currDataSize = size * size;
            } else {
                // isPacked => is symmetric
                currDataSize = ( ( 1.0 + ( double ) size ) *
                        ( ( double ) size ) / 2.0 );
            }
            long dataInBytes = (currDataSize + dataSize) * sizeof(double);
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

        this->DirichletPacks[i].Resize( matrixPerPack[i], dataSize );
        this->DirichletPacks[i].setMICratio( MICr );

        if ( configuration.load_balancing ) {
            this->DirichletPacks[i].enableLoadBalancing();
        } else {
            this->DirichletPacks[i].disableLoadBalancing();
        }

        for ( eslocal j = 0; j < matrixPerPack[i]; ++j ) {
            this->DirichletPacks[ i ].PreparePack( j, domains[accPreconditioners[i].at(j)].B1t_Dir_perm_vec.size(),
                    domains[accPreconditioners[i].at(j)].B1t_Dir_perm_vec.size(), symmetric );
        }
    }


    for (eslocal mic = 0 ; mic < configuration.N_MICS ; ++mic ) {
        this->DirichletPacks[mic].AllocateVectors( );
//        this->DirichletPacks[mic].SetDevice( mic );
        this->DirichletPacks[mic].SetDevice( target );        

		#pragma omp parallel for
        for (eslocal j = 0; j < accPreconditioners[mic].size(); ++j ) {
            eslocal d = accPreconditioners[mic].at(j);
            SEQ_VECTOR <eslocal> perm_vec = domains[d].B1t_Dir_perm_vec;
            SEQ_VECTOR <eslocal> perm_vec_full ( instance->K[d].rows );
            SEQ_VECTOR <eslocal> perm_vec_diff ( instance->K[d].rows );

            SEQ_VECTOR <eslocal> I_row_indices_p (instance->K[d].nnz);
            SEQ_VECTOR <eslocal> J_col_indices_p (instance->K[d].nnz);

            for (eslocal i = 0; i < perm_vec.size(); i++) {
                perm_vec[i] = perm_vec[i] - 1;
            }

            for (eslocal i = 0; i < perm_vec_full.size(); i++) {
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

            SEQ_VECTOR <SEQ_VECTOR<eslocal >> vec_I1_i2(K_modif.rows, SEQ_VECTOR<eslocal >(2, 1));
            eslocal offset = K_modif.CSR_I_row_indices[0] ? 1 : 0;

            for (eslocal i = 0; i < K_modif.rows;i++){
                vec_I1_i2[i][0] = perm_vec_full[i];
                vec_I1_i2[i][1] = i; // position to create reverse permutation
            }

            std::sort(vec_I1_i2.begin(), vec_I1_i2.end(), [](const SEQ_VECTOR <eslocal >& a, const SEQ_VECTOR<eslocal>& b) { return a[0] < b[0]; });

            // permutations made on matrix in COO format
            K_modif.ConvertToCOO(0);
            eslocal I_index,J_index;
            bool unsymmetric=!SYMMETRIC_SYSTEM;
            for (eslocal i = 0;i<K_modif.nnz;i++){
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
            for (eslocal i = 0; i<K_modif.nnz;i++){
                K_modif.I_row_indices[i] = I_row_indices_p[i];
                K_modif.J_col_indices[i] = J_col_indices_p[i];
            }
            K_modif.ConvertToCSRwithSort(1);
            {
                if (environment->print_matrices) {
                    std::ofstream osS(Logging::prepareFile(d, "K_modif"));
                    osS << K_modif;
                    osS.close();
                }
            }


            // ------------------------------------------------------------------------------------------------------------------
            bool diagonalized_K_rr = configuration.preconditioner == ESPRESO_PRECONDITIONER::SUPER_DIRICHLET;
            //        PRECONDITIONER==NONE              - 0
            //        PRECONDITIONER==LUMPED            - 1
            //        PRECONDITIONER==WEIGHT_FUNCTION   - 2
            //        PRECONDITIONER==DIRICHLET         - 3
            //        PRECONDITIONER==SUPER_DIRICHLET   - 4
            //        
            //        When next line is uncomment, var. PRECONDITIONER==DIRICHLET and PRECONDITIONER==SUPER_DIRICHLET provide identical preconditioner.
            //        bool diagonalized_K_rr = false
            // ------------------------------------------------------------------------------------------------------------------

            eslocal sc_size = perm_vec.size();

            if (sc_size == instance->K[d].rows) {
                domains[d].Prec = instance->K[d];
                domains[d].Prec.ConvertCSRToDense(1);
                // if physics.K[d] does not contain inner DOF
            } else {

                if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET) {
                    SparseSolverCPU createSchur;
                    //          createSchur.msglvl=1;
                    eslocal sc_size = perm_vec.size();
                    createSchur.ImportMatrix_wo_Copy(K_modif);
                    createSchur.Create_SC(domains[d].Prec, sc_size,false);
                    domains[d].Prec.ConvertCSRToDense(1);
                }
                else
                {
                    SparseMatrix K_rr;
                    SparseMatrix K_rs;
                    SparseMatrix K_sr;
                    SparseMatrix KsrInvKrrKrs; 

                    eslocal i_start = 0;
                    eslocal nonsing_size = K_modif.rows - sc_size - i_start;
                    eslocal j_start = nonsing_size;

                    K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, nonsing_size,j_start,sc_size);

                    if (SYMMETRIC_SYSTEM){
                        K_rs.MatTranspose(K_sr);
                    }
                    else
                    {
                        K_sr.getSubBlockmatrix_rs(K_modif,K_sr,j_start,sc_size,i_start, nonsing_size);
                    }

                    domains[d].Prec.getSubDiagBlockmatrix(K_modif,domains[d].Prec,nonsing_size,sc_size);
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
                        for (eslocal i = 0; i < K_rs.rows; i++) {
                            for (eslocal j = K_rs.CSR_I_row_indices[i]; j < K_rs.CSR_I_row_indices[i + 1]; j++) {
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
                    domains[d].Prec.MatAddInPlace(KsrInvKrrKrs,'N',-1);
                    //          if (!diagonalized_K_rr){
                    //				    domains[d].Prec.ConvertCSRToDense(1);
                    //          }
                }

            }

            if (environment->print_matrices) {
                std::ofstream osS(Logging::prepareFile(d, "S"));
                SparseMatrix SC =  domains[d].Prec;
                if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET){
                    SC.ConvertDenseToCSR(1);
                }
                osS << SC;
                osS.close();
            }

            // insert matrix to MatrixPack and delete it
            double * matrixPointer = this->DirichletPacks[mic].getMatrixPointer( j );
            memcpy( matrixPointer, &( domains[ d ].Prec.dense_values[ 0 ] ), 
                    this->DirichletPacks[ mic ].getDataLength( j ) * sizeof( double ) );
            SEQ_VECTOR<double>().swap( domains[d].Prec.dense_values );
            ESINFO(PROGRESS3) << Info::plain() << ".";
        }
    }

#pragma omp parallel num_threads(configuration.N_MICS)
    {
        this->DirichletPacks[omp_get_thread_num()].CopyToMIC();
    }


    // finish the blocks held on CPU
	#pragma omp parallel for
    for (size_t j = 0; j < hostPreconditioners.size(); ++j ) {
        eslocal d = hostPreconditioners.at(j);
        SEQ_VECTOR <eslocal> perm_vec = domains[d].B1t_Dir_perm_vec;
        SEQ_VECTOR <eslocal> perm_vec_full ( instance->K[d].rows );
        SEQ_VECTOR <eslocal> perm_vec_diff ( instance->K[d].rows );

        SEQ_VECTOR <eslocal> I_row_indices_p (instance->K[d].nnz);
        SEQ_VECTOR <eslocal> J_col_indices_p (instance->K[d].nnz);

        for (eslocal i = 0; i < perm_vec.size(); i++) {
            perm_vec[i] = perm_vec[i] - 1;
        }

        for (eslocal i = 0; i < perm_vec_full.size(); i++) {
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

        SEQ_VECTOR <SEQ_VECTOR<eslocal >> vec_I1_i2(K_modif.rows, SEQ_VECTOR<eslocal >(2, 1));
        eslocal offset = K_modif.CSR_I_row_indices[0] ? 1 : 0;

        for (eslocal i = 0; i < K_modif.rows;i++){
            vec_I1_i2[i][0] = perm_vec_full[i];
            vec_I1_i2[i][1] = i; // position to create reverse permutation
        }

        std::sort(vec_I1_i2.begin(), vec_I1_i2.end(), [](const SEQ_VECTOR <eslocal >& a, const SEQ_VECTOR<eslocal>& b) { return a[0] < b[0]; });

        // permutations made on matrix in COO format
        K_modif.ConvertToCOO(0);
        eslocal I_index,J_index;
        bool unsymmetric=!SYMMETRIC_SYSTEM;
        for (eslocal i = 0;i<K_modif.nnz;i++){
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
        for (eslocal i = 0; i<K_modif.nnz;i++){
            K_modif.I_row_indices[i] = I_row_indices_p[i];
            K_modif.J_col_indices[i] = J_col_indices_p[i];
        }
        K_modif.ConvertToCSRwithSort(1);
        {
            if (environment->print_matrices) {
                std::ofstream osS(Logging::prepareFile(d, "K_modif"));
                osS << K_modif;
                osS.close();
            }
        }


        // ------------------------------------------------------------------------------------------------------------------
        bool diagonalized_K_rr = configuration.preconditioner == ESPRESO_PRECONDITIONER::SUPER_DIRICHLET;
        //        PRECONDITIONER==NONE              - 0
        //        PRECONDITIONER==LUMPED            - 1
        //        PRECONDITIONER==WEIGHT_FUNCTION   - 2
        //        PRECONDITIONER==DIRICHLET         - 3
        //        PRECONDITIONER==SUPER_DIRICHLET   - 4
        //        
        //        When next line is uncomment, var. PRECONDITIONER==DIRICHLET and PRECONDITIONER==SUPER_DIRICHLET provide identical preconditioner.
        //        bool diagonalized_K_rr = false
        // ------------------------------------------------------------------------------------------------------------------

        eslocal sc_size = perm_vec.size();

        if (sc_size == instance->K[d].rows) {
            domains[d].Prec = instance->K[d];
            domains[d].Prec.ConvertCSRToDense(1);
            // if physics.K[d] does not contain inner DOF
        } else {

            if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET) {
                SparseSolverCPU createSchur;
                //          createSchur.msglvl=1;
                eslocal sc_size = perm_vec.size();
                createSchur.ImportMatrix_wo_Copy(K_modif);
                createSchur.Create_SC(domains[d].Prec, sc_size,false);
                domains[d].Prec.ConvertCSRToDense(1);
            }
            else
            {
                SparseMatrix K_rr;
                SparseMatrix K_rs;
                SparseMatrix K_sr;
                SparseMatrix KsrInvKrrKrs; 

                eslocal i_start = 0;
                eslocal nonsing_size = K_modif.rows - sc_size - i_start;
                eslocal j_start = nonsing_size;

                K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, nonsing_size,j_start,sc_size);

                if (SYMMETRIC_SYSTEM){
                    K_rs.MatTranspose(K_sr);
                }
                else
                {
                    K_sr.getSubBlockmatrix_rs(K_modif,K_sr,j_start,sc_size,i_start, nonsing_size);
                }

                domains[d].Prec.getSubDiagBlockmatrix(K_modif,domains[d].Prec,nonsing_size,sc_size);
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
                    for (eslocal i = 0; i < K_rs.rows; i++) {
                        for (eslocal j = K_rs.CSR_I_row_indices[i]; j < K_rs.CSR_I_row_indices[i + 1]; j++) {
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
                domains[d].Prec.MatAddInPlace(KsrInvKrrKrs,'N',-1);
                //          if (!diagonalized_K_rr){
                //				    domains[d].Prec.ConvertCSRToDense(1);
                //          }
            }

        }

        if (environment->print_matrices) {
            std::ofstream osS(Logging::prepareFile(d, "S"));
            SparseMatrix SC =  domains[d].Prec;
            if (configuration.preconditioner == ESPRESO_PRECONDITIONER::DIRICHLET){
                SC.ConvertDenseToCSR(1);
            }
            osS << SC;
            osS.close();
        }

        ESINFO(PROGRESS3) << Info::plain() << ".";
    }

    ESINFO(PROGRESS3);   
}
