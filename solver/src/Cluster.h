#ifdef WIN32
	#include "stdafx.h"
#endif

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
#include "FEM_Assembler.h"
#include "SparseSolver.h"
#include "TimeEval.h"
#include "Domain.h"
#include "DenseMatrixPack.h"

#include "utils.h"

#pragma once

class Cluster
{

public:
	// Constructor
	Cluster(int cluster_index);
	Cluster();

	// Cluster specific variables
	int cluster_global_index;
	int USE_DYNAMIC;
	int USE_KINV;
	int USE_HFETI;
	int SUBDOM_PER_CLUSTER;
	int NUMBER_OF_CLUSTERS;
	int DOFS_PER_NODE;

	int dual_size;
	string data_directory;

	// List of Domains
	SEQ_VECTOR <int>	domains_in_global_index;
	PAR_VECTOR <Domain> domains;

	int x_clust_size;
	SEQ_VECTOR <int> x_clust_domain_map_vec;

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
	int NUM_MICS;

	SparseSolver F0;
	SparseSolver F0_fast;
	SparseSolver Sa;

	SEQ_VECTOR <double> vec_d;

	SEQ_VECTOR <double> vec_b;
	SEQ_VECTOR <double> vec_b_compressed;

	SEQ_VECTOR <double>	vec_lambda;
	SEQ_VECTOR <double>	vec_alfa;
	SEQ_VECTOR <double>	vec_g0;
	SEQ_VECTOR <double>	vec_e0;

	SEQ_VECTOR<SEQ_VECTOR<double> > tm1;
	SEQ_VECTOR<SEQ_VECTOR<double> > tm2;
	SEQ_VECTOR<SEQ_VECTOR<double> > tm3;

	TimeEval cluster_time;

	// Functions of the class
	void LoadCluster(string directory_path, int use_dynamic_1_no_dynamic_0, int use_kinv_1_no_kinv_0 );

	void InitClusterPC ( int * subdomains_global_indices, int number_of_subdomains );
	void SetClusterPC  ( SEQ_VECTOR <SEQ_VECTOR <int> > & lambda_map); //, SEQ_VECTOR < int > & neigh_domains  );
	void SetClusterPC_AfterKplus ();
	void SetClusterHFETI ();

	void multKplusGlobal     ( SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, SEQ_VECTOR<int> & cluster_map_vec);
	void multKplusGlobal_l   ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in ); //, vector <double> & y_out, vector<int> & cluster_map_vec);
	void multKplusGlobal_Kinv( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );

	void CompressB0();
	void CreateG0();
	void CreateF0();
	void CreateSa();

	void Create_G1_perCluster();
	void CreateVec_d_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f );
	void CreateVec_b_perCluster( SEQ_VECTOR<SEQ_VECTOR <double> > & f );

	void Create_Kinv_perDomain();
	void Create_SC_perDomain();

	void B1_comp_MatVecSum(SEQ_VECTOR < SEQ_VECTOR <double> > & x_in, SEQ_VECTOR <double> & y_out, char T_for_transpose_N_for_non_transpose );

	void ShowTiming();
	void compress_lambda_vector(SEQ_VECTOR <double> & decompressed_vec_lambda);
	void decompress_lambda_vector(SEQ_VECTOR <double> & compressed_vec_lambda);

	//void GetProcessMemoryStat ( );
	//void GetMemoryStat( );

	SEQ_VECTOR < SEQ_VECTOR <int> >		my_comm_lambdas_indices;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_comm_lambdas;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_recv_lambdas;

	SEQ_VECTOR < SEQ_VECTOR <int> >		my_comm_lambdas_indices_comp;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_comm_lambdas_comp;
	SEQ_VECTOR < SEQ_VECTOR <double> >  my_recv_lambdas_comp;

	SEQ_VECTOR <int> my_neighs;
	SEQ_VECTOR <int> my_lamdas_indices;
	map <int,int> my_lamdas_map_indices;

	SEQ_VECTOR <int> my_lamdas_ddot_filter;
	SEQ_VECTOR <int> lambdas_filter;

	SEQ_VECTOR <double> compressed_tmp;
	SEQ_VECTOR <double> compressed_tmp2;

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

	int iter_cnt_comm;

	// variables for dynamic
	double dynamic_timestep;
	double dynamic_beta;
	double dynamic_gama;

	// dynamic
	void SetDynamicParameters(double set_dynamic_timestep, double set_dynamic_beta, double set_dynamic_gama);
};
