
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

#include "SparseMatrix.h"
#include "SparseSolver.h"
#include "Domain.h"
#include "DenseMatrixPack.h"

#include "utils.h"
#include "esconfig.h"
#include "esbasis.h"

#pragma once

class Cluster
{

public:
	// Constructor
	Cluster(eslocal cluster_index);
	Cluster();

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

	SparseSolver F0;
	SparseSolver F0_fast;
	SparseSolver Sa;
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
	void SetClusterPC  ( SEQ_VECTOR <SEQ_VECTOR <eslocal> > & lambda_map); //, SEQ_VECTOR < eslocal > & neigh_domains  );
	void SetClusterPC_AfterKplus ();
	void SetClusterHFETI ( bool R_from_mesh );

	void multKplusGlobal     ( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, SEQ_VECTOR<eslocal> & cluster_map_vec);
	void multKplusGlobal_l   ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in ); //, vector <double> & y_out, vector<eslocal> & cluster_map_vec);

	void multKplusGlobal_Kinv  ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );
	void multKplusGlobal_Kinv_2( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );

	void CompressB0();
	void CreateG0();
	void CreateF0();
	void CreateSa();

	void Create_G1_perCluster();
	void CreateVec_d_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f );
	void CreateVec_b_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f );

	void Create_Kinv_perDomain();
	void Create_SC_perDomain( bool USE_FLOAT );

	void B1_comp_MatVecSum(SEQ_VECTOR < SEQ_VECTOR <double> > & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose );

	void ShowTiming();
	void compress_lambda_vector(SEQ_VECTOR <double> & decompressed_vec_lambda);
	void decompress_lambda_vector(SEQ_VECTOR <double> & compressed_vec_lambda);

	//void GetProcessMemoryStat ( );
	//void GetMemoryStat( );

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
