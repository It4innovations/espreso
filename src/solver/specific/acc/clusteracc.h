
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
	ClusterAcc(const ESPRESOSolver &configuration, Instance *instance_in): ClusterBase(configuration, instance_in)
	{
			this->deleteMatrices = false;
			this->NUM_MICS = configuration.N_MICS;
			SetAcceleratorAffinity();
			ESINFO( DETAILS ) << "MICS: " << this->NUM_MICS;
	}

    virtual ~ClusterAcc();

    // creates Schur complement matrices
	void Create_SC_perDomain( bool USE_FLOAT );

    // factorizes stiffness matrices
	void SetupKsolvers ( );

    // assembles Dirichlet preconditioner
    void CreateDirichletPrec( Instance *instance );

    void multKplusGlobal_l_Acc(SEQ_VECTOR<SEQ_VECTOR<double> > & x_in, 
        double & CPUtime, 
        double * MICtime );

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

    eslocal NUM_MICS;

    bool deleteMatrices;
};

}


#endif /* SOLVER_SPECIFIC_ACC_CLUSTERACC_H_ */
