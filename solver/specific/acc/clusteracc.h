
#ifndef SOLVER_SPECIFIC_ACC_CLUSTERACC_H_
#define SOLVER_SPECIFIC_ACC_CLUSTERACC_H_

#include "../cluster.h"

#pragma offload_attribute(push,target(mic))
#include <unistd.h>
#pragma offload_attribute(pop))

class ClusterAcc: public ClusterBase
{

public:
	// Constructor
	ClusterAcc(eslocal cluster_index): ClusterBase(cluster_index) {
        this->deleteMatrices = false;
        this->NUM_MICS = 2;
    };
	ClusterAcc(): ClusterBase() {
        this->deleteMatrices = false;
        this->NUM_MICS = 2;
    };

    virtual ~ClusterAcc();

	void Create_SC_perDomain( bool USE_FLOAT );
    void Create_Kinv_perDomain();
	void SetupKsolvers ( );

//private:

    // packed matrices
    SEQ_VECTOR<DenseMatrixPack> B1KplusPacks;

    // number of accelerators
    eslocal NUM_MICS;

    // global solver for offloading all domains to Xeon Phi
    SEQ_VECTOR<SparseSolverAcc> solver;

    // array of matrix pointers per accelerator
    SEQ_VECTOR<SparseMatrix **> matricesPerAcc;

    // vector of length N_MIC of vectors of indices of subdomains on MICs
    SEQ_VECTOR<SEQ_VECTOR<eslocal> > accDomains;

    // vector of indices of domains on the host
    SEQ_VECTOR<eslocal> hostDomains;

    bool deleteMatrices;
};



#endif /* SOLVER_SPECIFIC_ACC_CLUSTERACC_H_ */
