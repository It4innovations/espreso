
#ifndef SOLVER_SPECIFIC_CLUSTER_H_
#define SOLVER_SPECIFIC_CLUSTER_H_


#include "feti/generic/SparseMatrix.h"
#include "feti/generic/Domain.h"
#include "densesolvers.h"

#include "feti/generic/timeeval.h"
#include "feti/generic/utils.h"
//#include "assembler/physics/assembler.h"

    namespace espreso {

//    struct Instance;


    class Domain;

    class ClusterBase
    {

    public:
        // Constructor
//        ClusterBase(const FETISolverConfiguration &configuration, esint cluster_index):
//        	configuration(configuration),
//            cluster_time("Cluster Timing "),
//            vec_fill_time("Reseting vec_g0 and vec_e0"),
//            loop_1_1_time("Loop 1: Kplus-sv, B0-mv, KpluR-mv"),
//            loop_1_2_time("Loop 1: vec_e0 and vec_g0"),
//            clusCP_time("Cluster CP - F0,GO,Sa,G0t,F0 "),
//            clus_F0_1_time("F0 solve - 1st "),
//            clus_F0_2_time("F0 solve - 2nd "),
//            clus_G0_time("G0  Mult "),
//            clus_G0t_time("G0t Mult "),
//            clus_Sa_time("Sa solve "),
//            loop_2_1_time("Loop2: Kplus-sv, B0-mv, Kplus-mv")
//        {
//            cluster_global_index = cluster_index;
//            iter_cnt_comm = 0;
//        };

//        ClusterBase(const FETISolverConfiguration &configuration):
//        	configuration(configuration),
//            cluster_time("Cluster Timing "),
//
//            vec_fill_time("Reseting vec_g0 and vec_e0"),
//            loop_1_1_time("Loop 1: Kplus-sv, B0-mv, KpluR-mv"),
//            loop_1_2_time("Loop 1: vec_e0 and vec_g0"),
//
//            clusCP_time("Cluster CP - F0,GO,Sa,G0t,F0 "),
//            clus_F0_1_time("F0 solve - 1st "),
//            clus_F0_2_time("F0 solve - 2nd "),
//            clus_G0_time("G0  Mult "),
//            clus_G0t_time("G0t Mult "),
//            clus_Sa_time("Sa solve "),
//
//            loop_2_1_time("Loop2: Kplus-sv, B0-mv, Kplus-mv")
//        {
//            iter_cnt_comm = 0;
//        }


        ClusterBase(const FETIConfiguration &configuration, DataHolder *instance_in):
        	configuration(configuration),
			instance(instance_in),

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
            ;
        }


        virtual ~ClusterBase() {};

        ClusterBase(const ClusterBase &other):
        	configuration(other.configuration),
        	instance(other.instance),

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

            loop_2_1_time	("Loop2: Kplus-sv, B0-mv, Kplus-mv"){

        }


        const FETIConfiguration &configuration;
        DataHolder *instance;

        // Cluster specific variables
        esint cluster_global_index;
        esint cluster_local_index;
        esint min_numClusters_per_MPI;


        esint USE_KINV;
        esint USE_HFETI;
        esint PAR_NUM_THREADS;
        esint SOLVER_NUM_THREADS;
        bool 	SYMMETRIC_SYSTEM;
        MatrixType mtype;

        esint dual_size;
        string data_directory;

        // List of Domains
        SEQ_VECTOR <esint> domains_in_global_index;
        PAR_VECTOR <Domain>  domains;

        esint x_clust_size;
        SEQ_VECTOR <esint> x_clust_domain_map_vec;

        // Global Matrices distributed per cluster
        SparseMatrix G1;
        SparseMatrix G1_comp;
        SparseMatrix G1t_comp;

        SparseMatrix G2;
        SparseMatrix G2_comp;
        SparseMatrix G2t_comp;


        SparseMatrix GGtinvM;
        SEQ_VECTOR <double> GGtinvV;

        // Matrices and solvers for Conj. projector
        SparseMatrix C;
        SparseMatrix Ct;
        DenseSolverMKL  CFCt;


        // Matrices and vectors of the cluster
        SparseMatrix G0;
        SparseMatrix G02;

        SparseMatrix F0_Mat;
        SparseMatrix B0Kplus;

        SparseSolverCPU F0;
        SparseSolverCPU F0_fast;

        SparseSolverCPU Sa;

        DenseSolverCPU 	Sa_dense_cpu;
        DenseSolverAcc  Sa_dense_acc;
        //SparseMatrix    SaMat;

        SEQ_VECTOR <double> vec_d;

        SEQ_VECTOR <double> vec_b;
        SEQ_VECTOR <double> vec_b_compressed;

        SEQ_VECTOR <double> vec_c;
        SEQ_VECTOR <double> vec_c_compressed;


        SEQ_VECTOR <double>	vec_lambda;
        SEQ_VECTOR <double>	vec_alfa;
        SEQ_VECTOR <double>	vec_g0;
        SEQ_VECTOR <double>	vec_e0;

        SEQ_VECTOR<SEQ_VECTOR<double> > tm1;
        SEQ_VECTOR<SEQ_VECTOR<double> > tm2;
        SEQ_VECTOR<SEQ_VECTOR<double> > tm3;

        TimeEval cluster_time;

        // Functions of the class

        void InitClusterPC   ( esint * subdomains_global_indices, esint number_of_subdomains );
        void SetClusterPC    (); //SEQ_VECTOR <SEQ_VECTOR <esint> > & lambda_map);
        void SetClusterHFETI ();

        virtual void SetupKsolvers ( ) = 0;
        void SetupPreconditioner ( );

        void multKplusGlobal     ( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, SEQ_VECTOR<esint> & cluster_map_vec);
        void multKplusGlobal_l   ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );

        void multKplusGlobal_Kinv  ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );
        void multKplusGlobal_Kinv_2( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );

        void CompressB0();
        void CreateG0();
        void CreateF0();
        void CreateSa();

//        void Create_G1_perCluster();
        void Create_G_perCluster();
        void Create_G_perSubdomain ( SparseMatrix &R_in, SparseMatrix &B_in, SparseMatrix &G_out );


        void Compress_G1();
        void Compress_G(SparseMatrix &G_in, SparseMatrix &G_comp_out);
        void CreateVec_d_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f );
        void CreateVec_b_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f );

        void CreateVec_c_perCluster( SEQ_VECTOR <double> & vec_c_out );
        void CreateVec_lb_perCluster( SEQ_VECTOR <double> & vec_lb_out );


        void Create_Kinv_perDomain();
        virtual void Create_SC_perDomain( bool USE_FLOAT ) = 0;
        virtual void CreateDirichletPrec( DataHolder *instance );

	void CreateDirichletPrec_perDomain( DataHolder *instance, size_t d);

	void B1_comp_MatVecSum(SEQ_VECTOR < SEQ_VECTOR <double> > & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose );

	void ShowTiming();
	void compress_lambda_vector(SEQ_VECTOR <double> & decompressed_vec_lambda);
	void decompress_lambda_vector(SEQ_VECTOR <double> & compressed_vec_lambda);

	SEQ_VECTOR < SEQ_VECTOR <esint> >		my_comm_lambdas_indices;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_comm_lambdas;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_recv_lambdas;

	SEQ_VECTOR < SEQ_VECTOR <esint> >		my_comm_lambdas_indices_comp;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_comm_lambdas_comp;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_recv_lambdas_comp;

	SEQ_VECTOR <esint> my_neighs;
	SEQ_VECTOR <esint> my_lamdas_indices;

	map <esint,esint> _my_lamdas_map_indices;

	SEQ_VECTOR <esint> my_lamdas_ddot_filter;
	SEQ_VECTOR <esint> lambdas_filter;

	SEQ_VECTOR <double> compressed_tmp;
	//SEQ_VECTOR <double> compressed_tmp2;

	SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster1;
	SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster2;
	SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster3;

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

};

}

#endif /* SOLVER_SPECIFIC_CLUSTER_H_ */
