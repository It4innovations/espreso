
#ifndef SOLVER_SPECIFIC_CLUSTER_H_
#define SOLVER_SPECIFIC_CLUSTER_H_


#include <omp.h>
#include "mpi.h"
#include "mkl.h"

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <map>

using std::vector;
using std::cout;
using std::map;
using std::make_pair;

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "../generic/SparseMatrix.h"
#include "../generic/Domain.h"
#include "acc/DenseMatrixPack.h"

#include "../generic/utils.h"
#include "esconfig.h"
#include "esbasis.h"

class Domain;


class ClusterBase
{

public:
	// Constructor
	ClusterBase(eslocal cluster_index):
		cluster_time("Cluster Timing "),
		vec_fill_time("Reseting vec_g0 and vec_e0"),
		loop_1_1_time("Loop 1: Kplus-sv, B0-mv, KpluR-mv"),
		loop_1_2_time("Loop 1: vec_e0 and vec_g0"),
		clusCP_time("Cluster CP - F0,GO,Sa,G0t,F0 "),
		clus_F0_1_time("F0 solve - 1st "),
		clus_F0_2_time("F0 solve - 2nd "),
		clus_G0_time("G0  Mult "),
		clus_G0t_time("G0t Mult "),
		clus_Sa_time("Sa solve "),
		loop_2_1_time("Loop2: Kplus-sv, B0-mv, Kplus-mv")
	{
		cluster_global_index = cluster_index;
		iter_cnt_comm = 0;
	};

	ClusterBase():
		cluster_time("Cluster Timing "),

		vec_fill_time("Reseting vec_g0 and vec_e0"),
		loop_1_1_time("Loop 1: Kplus-sv, B0-mv, KpluR-mv"),
		loop_1_2_time("Loop 1: vec_e0 and vec_g0"),

		clusCP_time("Cluster CP - F0,GO,Sa,G0t,F0 "),
		clus_F0_1_time("F0 solve - 1st "),
		clus_F0_2_time("F0 solve - 2nd "),
		clus_G0_time("G0  Mult "),
		clus_G0t_time("G0t Mult "),
		clus_Sa_time("Sa solve "),

		loop_2_1_time("Loop2: Kplus-sv, B0-mv, Kplus-mv")
	{
		iter_cnt_comm = 0;
	}

	virtual ~ClusterBase() {};

	// Cluster specific variables
	eslocal cluster_global_index;
	eslocal USE_DYNAMIC;
	eslocal USE_KINV;
	eslocal USE_HFETI;
	eslocal SUBDOM_PER_CLUSTER;
	eslocal NUMBER_OF_CLUSTERS;
	eslocal DOFS_PER_NODE;
	eslocal PAR_NUM_THREADS;
	eslocal SOLVER_NUM_THREADS;

	eslocal dual_size;
	string data_directory;

	// List of Domains
	SEQ_VECTOR <eslocal>	domains_in_global_index;
	PAR_VECTOR <Domain> domains;

	eslocal x_clust_size;
	SEQ_VECTOR <eslocal> x_clust_domain_map_vec;

	// Global Matrices distributed per cluster
	SparseMatrix G1;
	SparseMatrix G1_comp;
	SparseMatrix G1t_comp;
	SparseMatrix GGtinvM;
	SEQ_VECTOR <double> GGtinvV;

	// Matrices and vectors of the cluster
	SparseMatrix G0;
	SparseMatrix F0_Mat;
	SparseMatrix B0Kplus;

	// Packed matrices (mainly for MIC computation)
	SEQ_VECTOR <DenseMatrixPack> B1KplusPacks;

	// number of MIC
	eslocal NUM_MICS;

	SparseSolverCPU F0;
	SparseSolverCPU F0_fast;
	SparseSolverCPU Sa;
	SparseMatrix SaMat;

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

	void InitClusterPC ( eslocal * subdomains_global_indices, eslocal number_of_subdomains );
	void SetClusterPC  ( SEQ_VECTOR <SEQ_VECTOR <eslocal> > & lambda_map);
	void SetClusterPC_AfterKplus ();
	void SetClusterHFETI ( bool R_from_mesh );

	virtual void SetupKsolvers ( ) = 0;
	void ImportKmatrixAndRegularize ( SEQ_VECTOR <SparseMatrix> & K_in, const SEQ_VECTOR < SEQ_VECTOR < eslocal >> & fix_nodes );

	void multKplusGlobal     ( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, SEQ_VECTOR<eslocal> & cluster_map_vec);
	void multKplusGlobal_l   ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );

	void multKplusGlobal_Kinv  ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );
	void multKplusGlobal_Kinv_2( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );

	void CompressB0();
	void CreateG0();
	void CreateF0();
	void CreateSa();

	void Create_G1_perCluster();
	void Compress_G1();
	void CreateVec_d_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f );
	void CreateVec_b_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f );

	void Create_Kinv_perDomain();
	virtual void Create_SC_perDomain( bool USE_FLOAT ) = 0;

	void B1_comp_MatVecSum(SEQ_VECTOR < SEQ_VECTOR <double> > & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose );

	void ShowTiming();
	void compress_lambda_vector(SEQ_VECTOR <double> & decompressed_vec_lambda);
	void decompress_lambda_vector(SEQ_VECTOR <double> & compressed_vec_lambda);

	SEQ_VECTOR < SEQ_VECTOR <eslocal> >		my_comm_lambdas_indices;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_comm_lambdas;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_recv_lambdas;

	SEQ_VECTOR < SEQ_VECTOR <eslocal> >		my_comm_lambdas_indices_comp;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_comm_lambdas_comp;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_recv_lambdas_comp;

	SEQ_VECTOR <eslocal> my_neighs;
	SEQ_VECTOR <eslocal> my_lamdas_indices;

	map <eslocal,eslocal> _my_lamdas_map_indices;

	SEQ_VECTOR <eslocal> my_lamdas_ddot_filter;
	SEQ_VECTOR <eslocal> lambdas_filter;

	SEQ_VECTOR <double> compressed_tmp;
	//SEQ_VECTOR <double> compressed_tmp2;

	SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster1;
	SEQ_VECTOR<SEQ_VECTOR<double> > x_prim_cluster2;

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

	eslocal iter_cnt_comm;

	// variables for dynamic
	double dynamic_timestep;
	double dynamic_beta;
	double dynamic_gama;

	// dynamic
	void SetDynamicParameters(double set_dynamic_timestep, double set_dynamic_beta, double set_dynamic_gama);
};




#endif /* SOLVER_SPECIFIC_CLUSTER_H_ */
