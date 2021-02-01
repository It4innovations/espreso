
#ifndef SRC_SOLVER_SPECIFIC_ACC_SUPERCLUSTERACC_H_
#define SRC_SOLVER_SPECIFIC_ACC_SUPERCLUSTERACC_H_

#include "supercluster.h"

namespace espreso {

class SuperClusterAcc : public SuperClusterBase 
{
    public: 

    SuperClusterAcc( const FETIConfiguration & configuration, DataHolder *instance_in ):
        SuperClusterBase( configuration, instance_in ), 
       		cluster_time	("Cluster Timing "),

        vec_fill_time	("Reseting vec_g0 and vec_e0"),
        loop_1_1_time	("Loop 1: Kplus-sv, B0-mv, KpluR-mv"),
        loop_1_2_time	("Loop 1: vec_e0 and vec_g0"),

        clusCP_time		("Cluster CP - F0,GO,Sa,G0t,F0 "),
        clus_F0_1_time	("F0 solve - 1st "),
        clus_F0_2_time	("F0 solve - 2nd "),
        clus_G0_time	("G0  Mult "),
        clus_G0t_time	("G0t Mult "),
        clus_Sa_time	("Sa solve "),

        loop_2_1_time	("Loop2: Kplus-sv, B0-mv, Kplus-mv")
    {
        if (instance_in != NULL) {
            init();
        }
	    this->deleteMatrices = false;
		this->NUM_MICS = configuration.n_mics;
		SetAcceleratorAffinity();
		//ESINFO( DETAILS ) << "Number of Xeon Phi coprocessors: " << this->NUM_MICS;
    }

    virtual ~SuperClusterAcc();
    
    void init();

    void Create_SC_perDomain( bool USE_FLOAT );

    void SetupKsolvers( );
    
    void SetupPreconditioner( );

    void CreateDirichletPrec( DataHolder * instance );

    void multKplusGlobal_l_Acc(SEQ_VECTOR<SEQ_VECTOR<double>* > & x_in, 
        double & CPUtime, 
        double * MICtime );

    void SetAcceleratorAffinity( );

    SEQ_VECTOR <double> vec_g0;
    SEQ_VECTOR <double> vec_e0;

    SEQ_VECTOR <SEQ_VECTOR <double>* > tm1;
    SEQ_VECTOR <SEQ_VECTOR <double>* > tm2;
    SEQ_VECTOR <SEQ_VECTOR <double>* > tm3;
    
    // Matrices and vectors of the cluster
    SparseMatrix G0;
    SparseMatrix G02;
    SparseSolverCPU F0;
    SparseSolverCPU Sa;
    
    DenseSolverCPU  Sa_dense_cpu;
    DenseSolverAcc  Sa_dense_acc;

    SEQ_VECTOR <double> vec_alfa;
    SEQ_VECTOR <double> vec_lambda;

    // variables for time measurements
	TimeEvent vec_fill_time;
	TimeEvent loop_1_1_time;
	TimeEvent loop_1_2_time;

	TimeEvent clusCP_time;
	TimeEvent clus_F0_1_time;
	TimeEvent clus_F0_2_time;
	TimeEvent clus_G0_time;
	TimeEvent clus_G0t_time;
	TimeEvent clus_Sa_time;

	TimeEvent loop_2_1_time;
    TimeEvent cluster_time;


    // packed matrices
    SEQ_VECTOR<DenseMatrixPack> B1KplusPacks;

    // packed matrices for sparse solve on MIC
    SEQ_VECTOR<SparseMatrixPack> SparseKPack;
    
    // global solver for offloading all domains to Xeon Phi
    SEQ_VECTOR<SparseSolverAcc> solver;

    // array of matrix pointers per accelerator
    SEQ_VECTOR<SparseMatrix **> matricesPerAcc;

    // vector of length N_MIC of vectors of indices of subdomains on MICs
    SEQ_VECTOR<SEQ_VECTOR<esint> > accDomains;

    // vector of indices of domains on the host
    SEQ_VECTOR<esint> hostDomains;

    // packed Dirichlet preconditioners
    SEQ_VECTOR<DenseMatrixPack> DirichletPacks;

    // vector of length N_MIC of vectors of indices of preconditioners on MICs
    SEQ_VECTOR<SEQ_VECTOR<esint> > accPreconditioners;

    // vector of indices of domains on the host
    SEQ_VECTOR<esint> hostPreconditioners;

    // number of MPI processes per node
    esint MPI_per_node;

    // number of MPI processes per accelerator
    esint MPI_per_acc;
    
    // number of accelerators per MPI process
    esint acc_per_MPI;

    // id of devices to offload to
    SEQ_VECTOR<esint> myTargets;

    // local rank on the accelerator
    esint acc_rank;

    esint NUM_MICS;

    bool deleteMatrices;

};

}

#endif /* SRC_SOLVER_SPECIFIC_ACC_SUPERCLUSTERACC_H_ */ 
