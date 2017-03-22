
#ifndef SOLVER_SPECIFIC_ACC_CLUSTERACC_H_
#define SOLVER_SPECIFIC_ACC_CLUSTERACC_H_

#include "../cluster.h"

#pragma offload_attribute(push,target(mic))
#include <unistd.h>
#pragma offload_attribute(pop))

namespace espreso {

class ClusterAcc: public ClusterBase
{

public:
	// Constructor
	ClusterAcc(eslocal cluster_index): ClusterBase(cluster_index) {
        this->deleteMatrices = false;
        this->SetAcceleratorAffinity();
    };
	ClusterAcc(): ClusterBase() {
        this->deleteMatrices = false;
        this->SetAcceleratorAffinity();
    };

    virtual ~ClusterAcc();

	void Create_SC_perDomain( bool USE_FLOAT );
    void Create_Kinv_perDomain();
	void SetupKsolvers ( );
    void CreateDirichletPrec( Physics &physics );

    // sets affinity of processes on accelerators
    void SetAcceleratorAffinity();

//private:

    // packed matrices
    SEQ_VECTOR<DenseMatrixPack> B1KplusPacks;

    // packed matrices for sparse solve on MIC
    SEQ_VECTOR<SparseMatrixPack> SparseKPack;
    
    // global solver for offloading all domains to Xeon Phi
    SEQ_VECTOR<SparseSolverAcc> solver;

    // array of matrix pointers per accelerator
    SEQ_VECTOR<SparseMatrix **> matricesPerAcc;

    // vector of length N_MIC of vectors of indices of subdomains on MICs
    SEQ_VECTOR<SEQ_VECTOR<eslocal> > accDomains;

    // vector of indices of domains on the host
    SEQ_VECTOR<eslocal> hostDomains;

    // packed Dirichlet preconditioners
    SEQ_VECTOR<DenseMatrixPack> DirichletPacks;

    // vector of length N_MIC of vectors of indices of preconditioners on MICs
    SEQ_VECTOR<SEQ_VECTOR<eslocal> > accPreconditioners;

    // vector of indices of domains on the host
    SEQ_VECTOR<eslocal> hostPreconditioners;

    // number of MPI processes per node
    eslocal MPI_per_node;

    // number of MPI processes per accelerator
    eslocal MPI_per_acc;
    
    // number of accelerators per MPI process
    eslocal acc_per_MPI;

    // id of devices to offload to
    SEQ_VECTOR<eslocal> myTargets;

    // local rank on the accelerator
    eslocal acc_rank;

    bool deleteMatrices;
};

}


#endif /* SOLVER_SPECIFIC_ACC_CLUSTERACC_H_ */
